#!/usr/bin/env python3
# ============================================================================
# ENHANCED narrowPeak motif finder with robust error handling & cool features
# ============================================================================
# Adapted from the BEDPE-pair motif finder, but built for single-region
# narrowPeak files (MACS2/HiChIP-peaks style), e.g.:
#   chr1  91017  91597  peak_1  1522  .  17.8709  155.787  152.259  233
#   columns: chrom start end name score strand signalValue pValue qValue summit
#
# Since a narrowPeak line is a single region (not an anchor pair), the
# "canonical pair" orientation logic (">..<") from the BEDPE version doesn't
# apply directly. Instead this script reports, per peak:
#   - number of motif hits (forward/reverse)
#   - best hit score, orientation, and position
#   - distance from best hit to the peak summit
#   - an orientation-bias probability (like side_strength, but forward vs
#     reverse within the single peak)
#   - in --prob mode: a most-probable-orientation call ("<", ">", or "."
#     for no clear preference near 0.5); peaks with no motif hits at all
#     (probability < 0, i.e. no CTCF) are dropped from the output entirely
#   - an orientation "type" string (--stats mode): ">", "<", "><" (both), or "None"
#   - an overall quality score
#
# Features:
#   - Comprehensive error handling & input validation
#   - Multi-motif support (scan for multiple motifs simultaneously)
#   - Motif quality metrics (score, orientation bias)
#   - Distance analysis (motif to summit distance)
#   - Batch processing for memory efficiency
#   - JSON output for downstream analysis
#   - Detailed statistics & summary report
#   - Logging with multiple verbosity levels

import argparse
import time
import json
import sys
import logging
from pathlib import Path
from Bio import SeqIO, motifs
from sortedcontainers import SortedList
import numpy as np
from tqdm import tqdm
from collections import defaultdict
from dataclasses import dataclass

# ---- Setup logging ----
def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)

logger = None

# ---- Data structures ----
@dataclass
class MotifHit:
    """Single motif hit on a sequence"""
    position: int
    orientation: int  # +1 or -1
    score: float
    motif_id: str = ""

# ---- Utility functions ----
def order_key(x):
    """Key for sorting peaks by genomic position"""
    if isinstance(x, dict):
        return (x['chrom'], x['start'], x['end'])
    return (x.chrom, x.start, x.end)

def file_exists(file_path, context=""):
    """Validate file existence"""
    if not Path(file_path).exists():
        raise FileNotFoundError(f"{context}: {file_path} not found")
    return True

def validate_narrowpeak_format(line, line_num):
    """Validate narrowPeak line format (minimum 6 columns, ideally 10)"""
    vals = line.strip().split('\t')
    if len(vals) < 6:
        raise ValueError(f"Line {line_num}: Expected >=6 narrowPeak fields, got {len(vals)}")
    try:
        int(vals[1]), int(vals[2])
    except ValueError as e:
        raise ValueError(f"Line {line_num}: Invalid coordinate format - {e}")
    return vals

def normalize_chrom(chrom):
    """Normalize chromosome names"""
    chrom = str(chrom).strip()
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return chrom

# ---- Loading functions ----
def load_peaks(file_path, max_length=0, min_score=0, min_qvalue=None):
    """Load narrowPeak file with validation"""
    logger.info(f"Loading peaks from {file_path}")
    file_exists(file_path, "narrowPeak file")

    peaks = SortedList(key=order_key)
    skipped = 0

    try:
        with open(file_path) as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith("#") or line.startswith("track") or not line.strip():
                    continue
                try:
                    vals = validate_narrowpeak_format(line, line_num)

                    chrom = normalize_chrom(vals[0])
                    start, end = int(vals[1]), int(vals[2])
                    name = vals[3] if len(vals) > 3 else "."
                    score = float(vals[4]) if len(vals) > 4 and vals[4] != "." else 0.0
                    strand = vals[5] if len(vals) > 5 else "."
                    signal_value = float(vals[6]) if len(vals) > 6 and vals[6] != "." else None
                    p_value = float(vals[7]) if len(vals) > 7 and vals[7] != "." else None
                    q_value = float(vals[8]) if len(vals) > 8 and vals[8] != "." else None
                    summit_rel = int(vals[9]) if len(vals) > 9 and vals[9] != "." else (end - start) // 2

                    # Validation checks
                    if max_length > 0 and (end - start) > max_length:
                        skipped += 1
                        continue
                    if score < min_score:
                        skipped += 1
                        continue
                    if min_qvalue is not None and q_value is not None and q_value < min_qvalue:
                        skipped += 1
                        continue

                    peak = {
                        'chrom': chrom, 'start': start, 'end': end,
                        'name': name, 'score': score, 'strand': strand,
                        'signal_value': signal_value, 'p_value': p_value,
                        'q_value': q_value, 'summit_rel': summit_rel,
                        'summit_abs': start + summit_rel
                    }
                    peaks.add(peak)
                except ValueError as e:
                    logger.warning(f"Skipping line {line_num}: {e}")
                    skipped += 1
    except IOError as e:
        raise IOError(f"Cannot read narrowPeak file: {e}")

    logger.info(f"Loaded {len(peaks)} peaks (skipped {skipped})")
    return peaks

def load_genome(fasta_path, batch_mode=False):
    """Load genome FASTA into dictionary (or return path for batch mode)"""
    file_exists(fasta_path, "Genome FASTA")

    if batch_mode:
        logger.info(f"Genome file: {fasta_path} (batch mode - will stream sequences)")
        return fasta_path

    logger.info(f"Loading genome from {fasta_path}")
    genome = {}
    filtered = 0

    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            if any(x in record.id for x in ["_", "*", "EBV", "decoy", "random"]):
                filtered += 1
                continue
            genome[record.id] = record.seq
        logger.info(f"Loaded {len(genome)} contigs (filtered {filtered})")
    except Exception as e:
        raise IOError(f"Error parsing FASTA: {e}")

    if not genome:
        raise ValueError("No valid contigs found in genome file")

    return genome

def get_sequence_batch(genome_path, chrom, start, end):
    """Get sequence from genome file (for batch mode)"""
    for record in SeqIO.parse(genome_path, "fasta"):
        if record.id == chrom:
            return record.seq[max(0, start):end]
    return None

def load_motif_pfm(motif_path):
    """Load motif PFM file safely"""
    file_exists(motif_path, "Motif PFM")
    logger.info(f"Loading motif from {motif_path}")

    try:
        with open(motif_path) as f:
            pfm = motifs.read(f, "pfm")
            pssm = pfm.pssm
            logger.info(f"Motif loaded: {len(pssm)} positions")
            return pssm
    except Exception as e:
        raise ValueError(f"Error loading motif: {e}")

def load_multiple_motifs(motif_paths):
    """Load multiple motif files"""
    logger.info(f"Loading {len(motif_paths)} motifs")
    motifs_dict = {}

    for i, path in enumerate(motif_paths):
        try:
            pssm = load_motif_pfm(path)
            motif_id = f"motif_{i+1}_{Path(path).stem}"
            motifs_dict[motif_id] = pssm
        except Exception as e:
            logger.error(f"Failed to load {path}: {e}")

    if not motifs_dict:
        raise ValueError("No motifs loaded successfully")

    return motifs_dict

# ---- Motif scanning ----
def scan_sequence_raw(pssm, seq, threshold=7.0):
    """Scan sequence and return raw pssm.search results (position, score)"""
    try:
        return list(pssm.search(seq, threshold=threshold))
    except Exception as e:
        logger.debug(f"Error scanning sequence: {e}")
        return []

def scan_sequence(pssm, seq, threshold=7.0):
    """Scan sequence and return all hits above threshold as MotifHit objects"""
    hits = []
    try:
        results = scan_sequence_raw(pssm, seq, threshold)
        for position, score in results:
            orientation = 1 if position >= 0 else -1
            hits.append(MotifHit(
                position=abs(position),
                orientation=orientation,
                score=score
            ))
    except Exception as e:
        logger.debug(f"Error scanning sequence: {e}")
    return hits

def orientation_bias(results):
    """
    Probability that motif hits within the peak are forward-oriented,
    weighted by score (analogous to side_strength in the paired version).
    Returns -1 if no motifs found.
    """
    if not results:
        return -1.0

    forward = sum(2**score for orient, score in results if orient >= 0)
    reverse = sum(2**score for orient, score in results if orient < 0)
    total = forward + reverse

    if total == 0:
        return -1.0
    return forward / total

def scan_motif_probabilistic(pssm, seq, threshold=2.0):
    """Probabilistic scoring mode: forward-orientation bias (0-1) or -1 if none"""
    results = scan_sequence_raw(pssm, seq, threshold=threshold)
    return orientation_bias(results)

def classify_orientation(prob, band=0.1):
    """
    Map a forward-orientation probability to a most-probable-orientation call:
    - "NaN" if prob is negative (no motif found, i.e. the -1 sentinel)
    - "<"   if prob is close to 0 (reverse-dominant)
    - ">"   if prob is close to 1 (forward-dominant)
    - "."   if prob is close to 0.5 (no preferential orientation)
    `band` sets how close to 0.5 counts as "no preference".
    """
    if prob is None or prob < 0:
        return "NaN"
    if abs(prob - 0.5) <= band:
        return "."
    return ">" if prob > 0.5 else "<"

def scan_motif_stats(pssm, seq, threshold=7.0):
    """
    Stats/orientation string mode for a single region:
    - ">"    : forward-oriented hit(s) only
    - "<"    : reverse-oriented hit(s) only
    - "><"   : both orientations present
    - "None" : no motifs found
    """
    results = scan_sequence_raw(pssm, seq, threshold=threshold)
    if not results:
        return "None"

    has_fwd = any(orient >= 0 for orient, score in results)
    has_rev = any(orient < 0 for orient, score in results)

    if has_fwd and has_rev:
        return "><"
    if has_fwd:
        return ">"
    return "<"

def scan_motif_binary(pssm, seq, threshold=7.0):
    """
    Binary mode: True if the peak contains at least one motif hit
    above threshold, using original-style heuristics.
    """
    results = scan_sequence_raw(pssm, seq, threshold=threshold)
    return len(results) > 0

def calculate_motif_distance(hits, summit_abs, region_start):
    """Distance from the best-scoring motif hit to the peak summit (bp)"""
    if not hits:
        return None
    best_hit = max(hits, key=lambda h: h.score)
    hit_abs_pos = region_start + best_hit.position
    return abs(hit_abs_pos - summit_abs)

def calculate_quality_score(hits, at_summit_bonus=0.0):
    """
    Overall quality score combining:
    - Mean motif score
    - Number of hits (log-scaled bonus)
    - Bonus for proximity to summit (passed in by caller)
    """
    if not hits:
        return 0.0

    avg_score = np.mean([h.score for h in hits])
    count_bonus = min(np.log1p(len(hits)), 2.0)
    score = avg_score + count_bonus + at_summit_bonus
    return min(score, 20.0)

# ---- Main processing ----
def process_peaks(peaks, genome, pssms, threshold=7.0,
                   distance_threshold=None, probabilistic=False, stats_mode=False):
    """
    Process all peaks and scan for motifs.

    Modes:
    - Default: hit counts, best score/orientation, distance to summit, quality
    - probabilistic=True: forward-orientation bias (-1 to 1) per peak
    - stats_mode=True: orientation strings (">", "<", "><", "None")
    """
    logger.info(f"Scanning {len(peaks)} peaks with {len(pssms)} motif(s)")
    if probabilistic:
        logger.info("Using probabilistic scoring mode")
    if stats_mode:
        logger.info("Using stats/orientation mode")

    results = SortedList(key=order_key)
    stats = defaultdict(int)

    for peak in tqdm(peaks, desc="Scanning", ncols=90, colour="cyan"):
        try:
            chrom, start, end = peak['chrom'], peak['start'], peak['end']
            summit_abs = peak['summit_abs']

            # Extract sequence
            if isinstance(genome, str):  # Batch mode
                seq = get_sequence_batch(genome, chrom, start, end)
            else:
                seq = genome.get(chrom, None)

            if seq is None:
                logger.debug(f"Chromosome not found: {chrom}")
                stats['missing_chrom'] += 1
                continue

            seq = seq[start:end]
            if len(seq) == 0:
                stats['empty_seq'] += 1
                continue

            seq_str = str(seq)

            if probabilistic:
                for motif_id, pssm in pssms.items():
                    prob = scan_motif_probabilistic(pssm, seq_str, threshold=2.0)

                    # Negative prob means no motif hits at all (no CTCF present) -
                    # skip the peak entirely rather than writing a NaN row.
                    if prob < 0:
                        stats['no_motif'] += 1
                        continue

                    result = {
                        'chrom': chrom, 'start': start, 'end': end,
                        'name': peak['name'], 'score': peak['score'],
                        'prob_forward': prob,  # always 0-1 here; no-motif peaks are skipped above
                        'orientation': classify_orientation(prob),  # >, <, or .
                        'motif_id': motif_id
                    }
                    results.add(result)
                    stats['probabilistic'] += 1

            elif stats_mode:
                for motif_id, pssm in pssms.items():
                    orientation_type = scan_motif_stats(pssm, seq_str, threshold)
                    result = {
                        'chrom': chrom, 'start': start, 'end': end,
                        'name': peak['name'], 'score': peak['score'],
                        'type': orientation_type,
                        'motif_id': motif_id
                    }
                    results.add(result)
                    stats['stats'] += 1

            else:
                all_hits = []
                passed_filter = False

                for motif_id, pssm in pssms.items():
                    has_motif = scan_motif_binary(pssm, seq_str, threshold)
                    if has_motif:
                        passed_filter = True
                        hits = scan_sequence(pssm, seq_str, threshold)
                        for h in hits:
                            h.motif_id = motif_id
                        all_hits.extend(hits)

                if not passed_filter:
                    stats['no_motif'] += 1
                    continue

                dist = calculate_motif_distance(all_hits, summit_abs, start)

                if distance_threshold and dist is not None and dist > distance_threshold:
                    stats['far_from_summit'] += 1
                    continue

                at_summit_bonus = 1.0 if (dist is not None and dist <= 20) else 0.0
                quality = calculate_quality_score(all_hits, at_summit_bonus)

                best_hit = max(all_hits, key=lambda h: h.score)
                n_fwd = sum(1 for h in all_hits if h.orientation >= 0)
                n_rev = sum(1 for h in all_hits if h.orientation < 0)

                result = {
                    'chrom': chrom, 'start': start, 'end': end,
                    'name': peak['name'], 'score': peak['score'],
                    'motif_hits': len(all_hits),
                    'motif_hits_fwd': n_fwd,
                    'motif_hits_rev': n_rev,
                    'best_score': best_hit.score,
                    'best_orientation': best_hit.orientation,
                    'dist_to_summit': dist,
                    'quality_score': quality
                }
                results.add(result)
                stats['passed'] += 1

        except Exception as e:
            logger.debug(f"Error processing peak {peak}: {e}")
            stats['error'] += 1

    return results, stats

# ---- Output functions ----
def save_narrowpeak(filename, results, include_scores=False):
    """
    Save results to a narrowPeak-like format.
    Handles three output layouts:
    - Default: chrom/start/end/name/score + motif count + quality + distance
    - Probabilistic: chrom/start/end/name/score + prob_forward + orientation (>, <, or .) + motif_id
      (peaks with no motif hits at all - prob < 0 - are excluded upstream, since that means no CTCF)
    - Stats: chrom/start/end/name/score + orientation type + motif_id
    """
    logger.info(f"Writing {len(results)} results to {filename}")

    try:
        with open(filename, "w") as f:
            for r in results:
                line = f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['name']}\t{r['score']}"

                if 'prob_forward' in r:
                    line += f"\t{r['prob_forward']:.2f}\t{r['orientation']}"
                    if 'motif_id' in r:
                        line += f"\t{r['motif_id']}"

                elif 'type' in r:
                    line += f"\t{r['type']}"
                    if 'motif_id' in r:
                        line += f"\t{r['motif_id']}"

                elif include_scores and 'motif_hits' in r:
                    line += (f"\t{r['motif_hits']}\t{r['motif_hits_fwd']}\t{r['motif_hits_rev']}"
                             f"\t{r['best_score']:.2f}\t{r['best_orientation']}"
                             f"\t{r.get('dist_to_summit', '.')}\t{r['quality_score']:.2f}")

                f.write(line + "\n")
    except IOError as e:
        raise IOError(f"Cannot write output file: {e}")

def save_json(filename, results):
    """Save results as JSON for downstream analysis"""
    logger.info(f"Writing JSON results to {filename}")

    try:
        with open(filename, "w") as f:
            json.dump([dict(r) for r in results], f, indent=2, default=str)
    except Exception as e:
        raise IOError(f"Cannot write JSON file: {e}")

def generate_report(results, stats, output_file=None):
    """Generate summary statistics report"""
    report = []
    report.append("\n" + "="*60)
    report.append("NARROWPEAK MOTIF SCANNING REPORT")
    report.append("="*60)
    report.append(f"\nTotal peaks passed filters: {stats['passed']}")
    report.append(f"Skipped - no motif found: {stats['no_motif']}")
    report.append(f"Skipped - motif too far from summit: {stats['far_from_summit']}")
    report.append(f"Skipped - missing chromosome: {stats['missing_chrom']}")
    report.append(f"Skipped - empty sequence: {stats['empty_seq']}")
    report.append(f"Errors encountered: {stats['error']}")

    if results and stats['passed'] > 0:
        avg_quality = np.mean([r['quality_score'] for r in results if 'quality_score' in r])
        report.append(f"\nAverage quality score: {avg_quality:.2f}")

        hits_dist = [r['motif_hits'] for r in results if 'motif_hits' in r]
        if hits_dist:
            report.append(f"\nMotif hits per peak: mean={np.mean(hits_dist):.2f}, max={max(hits_dist)}")

        valid_dists = [r['dist_to_summit'] for r in results
                        if r.get('dist_to_summit') is not None]
        if valid_dists:
            report.append(f"\nDistance to summit (bp):")
            report.append(f"  Mean: {np.mean(valid_dists):.0f}")
            report.append(f"  Median: {np.median(valid_dists):.0f}")
            report.append(f"  Max: {max(valid_dists):.0f}")

    report.append("\n" + "="*60 + "\n")
    report_text = "\n".join(report)

    logger.info(report_text)

    if output_file:
        try:
            with open(output_file, "w") as f:
                f.write(report_text)
        except IOError as e:
            logger.warning(f"Cannot write report: {e}")

    return report_text

# ---- Main entrypoint ----
def main():
    parser = argparse.ArgumentParser(
        description="Enhanced narrowPeak motif finder with multi-motif support and robust validation"
    )
    parser.add_argument("-i", "--input", required=True, help="Input narrowPeak file")
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA file")
    parser.add_argument("-m", "--motif", required=True, nargs='+',
                       help="One or more motif PFM files")
    parser.add_argument("-o", "--output", default="output_motifs.narrowPeak",
                       help="Output narrowPeak-like file")
    parser.add_argument("--json", help="Output JSON file with full results")
    parser.add_argument("--report", help="Output statistics report")
    parser.add_argument("--threshold", type=float, default=7.0,
                       help="PSSM score threshold (default=7.0)")
    parser.add_argument("--max-distance", type=int,
                       help="Max distance from motif to summit (bp)")
    parser.add_argument("--min-score", type=float, default=0,
                       help="Minimum narrowPeak score (col 5)")
    parser.add_argument("--min-qvalue", type=float,
                       help="Minimum -log10(qvalue) (col 9)")
    parser.add_argument("--include-scores", action="store_true",
                       help="Include motif scores in output")
    parser.add_argument("--prob", action="store_true",
                       help="Probabilistic scoring mode (forward-orientation bias, -1 to 1)")
    parser.add_argument("--stats", action="store_true",
                       help="Stats mode (orientation strings like >, <, ><, or None)")
    parser.add_argument("--batch", action="store_true",
                       help="Batch mode: stream genome (slower, less memory)")
    parser.add_argument("--verbose", action="store_true",
                       help="Verbose logging")

    args = parser.parse_args()

    global logger
    logger = setup_logging(args.verbose)

    try:
        motif_files = args.motif if isinstance(args.motif, list) else [args.motif]

        logger.info("=" * 60)
        peaks = load_peaks(args.input, min_score=args.min_score, min_qvalue=args.min_qvalue)
        genome = load_genome(args.genome, batch_mode=args.batch)
        pssms = load_multiple_motifs(motif_files)

        logger.info("=" * 60)
        start = time.time()
        results, stats = process_peaks(
            peaks, genome, pssms,
            threshold=args.threshold,
            distance_threshold=args.max_distance,
            probabilistic=args.prob,
            stats_mode=args.stats
        )
        elapsed = time.time() - start

        logger.info("=" * 60)
        save_narrowpeak(args.output, results, include_scores=args.include_scores)
        if args.json:
            save_json(args.json, results)

        report_file = args.report if args.report else None
        generate_report(results, stats, report_file)

        logger.info(f"Completed in {elapsed:.2f}s")
        logger.info(f"Results saved to {args.output}")

    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=args.verbose)
        sys.exit(1)

if __name__ == "__main__":
    main()