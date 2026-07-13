#!/usr/bin/env python3
# ============================================================================
# BED Motif Finder - Simplified motif scanning for single intervals
# Keeps original algorithmic logic: binary (default), probabilistic (--prob), stats (--stats)
# ============================================================================

import argparse
import time
import json
import sys
import logging
from pathlib import Path
from Bio import SeqIO, motifs
from sortedcontainers import SortedList
import numpy as np
import vcf
from tqdm import tqdm
from collections import defaultdict

# ---- Logging ----
def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)

logger = None

# ---- File I/O & Validation ----
def file_exists(file_path, context=""):
    if not Path(file_path).exists():
        raise FileNotFoundError(f"{context}: {file_path} not found")
    return True

def normalize_chrom(chrom):
    chrom = str(chrom).strip()
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return chrom

def validate_bed_format(line, line_num):
    """Validate BED line format (at least 3 fields)"""
    vals = line.strip().split('\t')
    if len(vals) < 3:
        raise ValueError(f"Line {line_num}: Expected ≥3 BED fields, got {len(vals)}")
    try:
        int(vals[1]), int(vals[2])
    except ValueError as e:
        raise ValueError(f"Line {line_num}: Invalid coordinate format - {e}")
    return vals

# ---- Loading ----
def load_bed(bed_path):
    """Load BED file into list of dicts"""
    logger.info(f"Loading BED file from {bed_path}")
    file_exists(bed_path, "BED file")
    
    intervals = []
    skipped = 0
    
    try:
        with open(bed_path) as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith("#") or not line.strip():
                    continue
                try:
                    vals = validate_bed_format(line, line_num)
                    
                    chrom = normalize_chrom(vals[0])
                    start = int(vals[1])
                    end = int(vals[2])
                    name = vals[3] if len(vals) > 3 else f"{chrom}:{start}-{end}"
                    score = vals[4] if len(vals) > 4 else "."
                    strand = vals[5] if len(vals) > 5 else "."
                    
                    interval = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'name': name,
                        'score': score,
                        'strand': strand
                    }
                    intervals.append(interval)
                except ValueError as e:
                    logger.warning(f"Skipping line {line_num}: {e}")
                    skipped += 1
    except IOError as e:
        raise IOError(f"Cannot read BED file: {e}")
    
    logger.info(f"Loaded {len(intervals)} intervals (skipped {skipped})")
    return intervals

def load_genome(fasta_path, batch_mode=False):
    """Load genome FASTA"""
    file_exists(fasta_path, "Genome FASTA")
    
    if batch_mode:
        logger.info(f"Genome file: {fasta_path} (batch mode)")
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
    """Get sequence from genome file (batch mode)"""
    for record in SeqIO.parse(genome_path, "fasta"):
        if record.id == chrom:
            return record.seq[max(0, start):end]
    return None

def load_motif_pfm(motif_path):
    """Load motif PFM file"""
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

def load_snps(vcf_path):
    """Load SNPs from VCF"""
    file_exists(vcf_path, "VCF file")
    logger.info(f"Loading SNPs from {vcf_path}")
    
    snps = SortedList()
    count = 0
    
    try:
        reader = vcf.Reader(filename=vcf_path)
        for record in reader:
            if not record.CHROM or record.is_indel or record.is_sv:
                continue
            if any(x in record.CHROM for x in ["_", "*", "EBV", "decoy"]):
                continue
            
            chr_name = normalize_chrom(record.CHROM)
            pos = record.POS
            
            try:
                if len(record.ALT) > 1:
                    alleles = [str(record.ALT[0]), str(record.ALT[1])]
                else:
                    alleles = [str(record.REF), str(record.ALT[0])]
                
                snps.add((chr_name, pos, alleles))
                count += 1
            except Exception as e:
                logger.debug(f"Skipping SNP at {record.CHROM}:{pos}: {e}")
                continue
    except Exception as e:
        raise IOError(f"Error reading VCF file: {e}")
    
    logger.info(f"Loaded {count} SNPs")
    return snps if count > 0 else None

def update_sequence(seq, snps, offset):
    """Apply SNP alleles to sequence"""
    seq_list = list(str(seq))
    for chr_name, snp_pos, alleles in snps:
        rel_pos = snp_pos - offset
        if 0 <= rel_pos < len(seq_list):
            seq_list[rel_pos] = alleles[0]
    return ''.join(seq_list)

# ---- Motif Scanning (Original Logic Restored) ----
def scan_sequence_raw(pssm, seq, threshold=7.0):
    """Raw PSSM search results (position, score)"""
    try:
        return list(pssm.search(seq, threshold=threshold))
    except Exception as e:
        logger.debug(f"Error scanning: {e}")
        return []

def side_strength(results):
    """
    ORIGINAL: Calculate probabilistic orientation
    Returns -1 if no motifs, else 0-1 (left bias)
    """
    if not results:
        return -1.0
    
    left = sum(2**score for orient, score in results if orient < 0)
    right = sum(2**score for orient, score in results if orient >= 0)
    total = left + right
    
    if total == 0:
        return -1.0
    return left / total

def scan_motif_binary(pssm, seq, threshold=7.0):
    """
    ORIGINAL: Binary mode - returns True if motif found with good score
    """
    results = scan_sequence_raw(pssm, seq, threshold)
    return len(results) > 0

def scan_motif_probabilistic(pssm, seq, threshold=2.0):
    """
    ORIGINAL: Probabilistic mode - returns probability (-1 to 1)
    """
    results = scan_sequence_raw(pssm, seq, threshold)
    return side_strength(results)

def classify_orientation(prob, band=0.1):
    """
    Map a left-orientation probability to a most-probable-orientation call:
    - "NaN" if prob is negative (no motif found, i.e. the -1 sentinel)
    - "<"   if prob is close to 0 (right/forward-dominant)
    - ">"   if prob is close to 1 (left/reverse-dominant)
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
    ORIGINAL: Stats mode - returns orientation string
    ">" = forward orientation
    "<" = reverse orientation  
    "None" = no motif found
    """
    results = scan_sequence_raw(pssm, seq, threshold)
    
    if not results:
        return "None"
    
    # Get best hit orientation
    best_orient = max(results, key=lambda x: x[1])[0]
    return ">" if best_orient >= 0 else "<"

def count_motifs(pssm, seq, threshold=7.0):
    """Count number of motif hits"""
    return len(scan_sequence_raw(pssm, seq, threshold))

def get_best_score(pssm, seq, threshold=7.0):
    """Get best motif score"""
    results = scan_sequence_raw(pssm, seq, threshold)
    if not results:
        return None
    return max(results, key=lambda x: x[1])[1]

# ---- Main Processing ----
def process_bed(intervals, genome, pssm, snps=None, threshold=7.0,
                probabilistic=False, stats_mode=False):
    """
    Process BED intervals and scan for motifs
    
    Modes:
    - Default: Return only intervals with motifs (binary)
    - probabilistic=True: Return probability scores (0-1) + orientation call (>, <, .)
      Intervals with no motif hits at all (prob < 0) are skipped entirely.
    - stats_mode=True: Return orientation strings (">" / "<" / "None")
    """
    logger.info(f"Scanning {len(intervals)} intervals")
    if probabilistic:
        logger.info("Mode: Probabilistic scoring")
    elif stats_mode:
        logger.info("Mode: Orientation classification")
    else:
        logger.info("Mode: Binary filtering")
    
    results = []
    stats = defaultdict(int)
    
    for interval in tqdm(intervals, desc="Scanning", ncols=80, colour="cyan"):
        try:
            chrom = interval['chrom']
            start = interval['start']
            end = interval['end']
            
            # Get sequence
            if isinstance(genome, str):  # Batch mode
                seq = get_sequence_batch(genome, chrom, start, end)
            else:
                seq = genome.get(chrom, None)
            
            if seq is None:
                logger.debug(f"Chromosome not found: {chrom}")
                stats['missing_chrom'] += 1
                continue
            
            # Extract subsequence
            seq = seq[start:end]
            
            if len(seq) == 0:
                stats['empty_seq'] += 1
                continue
            
            # Apply SNPs if available
            if snps:
                snps_in_interval = [s for s in snps if s[0] == chrom and start <= s[1] <= end]
                if snps_in_interval:
                    seq = update_sequence(seq, snps_in_interval, start)
            
            seq_str = str(seq)
            
            # Process based on mode
            if probabilistic:
                # ORIGINAL: Probabilistic mode
                prob = scan_motif_probabilistic(pssm, seq_str, threshold=2.0)

                # Negative prob means no motif hits at all (no CTCF present) -
                # skip the interval entirely rather than writing a NaN row.
                if prob < 0:
                    stats['no_motif'] += 1
                    continue

                result = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': interval['name'],
                    'score': interval['score'],
                    'strand': interval['strand'],
                    'prob': prob,  # always 0-1 here; no-motif intervals are skipped above
                    'orientation': classify_orientation(prob)  # >, <, or .
                }
                results.append(result)
                stats['with_score'] += 1
                
            elif stats_mode:
                # ORIGINAL: Stats mode
                orientation = scan_motif_stats(pssm, seq_str, threshold)
                
                result = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': interval['name'],
                    'score': interval['score'],
                    'strand': interval['strand'],
                    'motif_orientation': orientation
                }
                results.append(result)
                stats['with_orientation'] += 1
                
            else:
                # ORIGINAL: Binary mode
                has_motif = scan_motif_binary(pssm, seq_str, threshold)
                
                if not has_motif:
                    stats['no_motif'] += 1
                    continue
                
                # Get extra metrics
                hit_count = count_motifs(pssm, seq_str, threshold)
                best_score = get_best_score(pssm, seq_str, threshold)
                
                result = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': interval['name'],
                    'score': interval['score'],
                    'strand': interval['strand'],
                    'motif_hits': hit_count,
                    'motif_score': best_score
                }
                results.append(result)
                stats['passed'] += 1
            
        except Exception as e:
            logger.debug(f"Error processing {interval['name']}: {e}")
            stats['error'] += 1
    
    return results, stats

# ---- Output ----
def save_bed(filename, results, mode='default'):
    """
    Save results to BED format
    
    Modes:
    - default: BED + motif_hits + motif_score
    - probabilistic: BED + probability score + orientation call (>, <, or .)
      (intervals with no motif hits at all - prob < 0 - are excluded upstream, since that means no CTCF)
    - stats: BED + orientation
    """
    logger.info(f"Writing {len(results)} results to {filename}")
    
    try:
        with open(filename, "w") as f:
            for r in results:
                line = f"{r['chrom']}\t{r['start']}\t{r['end']}"
                line += f"\t{r['name']}\t{r['score']}\t{r['strand']}"
                
                if mode == 'probabilistic':
                    line += f"\t{r['prob']:.2f}\t{r['orientation']}"
                elif mode == 'stats':
                    line += f"\t{r['motif_orientation']}"
                else:  # default
                    line += f"\t{r['motif_hits']}\t{r['motif_score']:.2f}"
                
                f.write(line + "\n")
    except IOError as e:
        raise IOError(f"Cannot write output: {e}")

def save_json(filename, results):
    """Save results as JSON"""
    logger.info(f"Writing JSON to {filename}")
    try:
        with open(filename, "w") as f:
            json.dump(results, f, indent=2, default=str)
    except Exception as e:
        raise IOError(f"Cannot write JSON: {e}")

def generate_report(results, stats, output_file=None):
    """Generate summary report"""
    report = []
    report.append("\n" + "="*60)
    report.append("BED MOTIF SCANNING REPORT")
    report.append("="*60)
    report.append(f"\nTotal intervals processed: {len(results)}")
    report.append(f"Skipped - no motif: {stats['no_motif']}")
    report.append(f"Skipped - missing chrom: {stats['missing_chrom']}")
    report.append(f"Skipped - empty sequence: {stats['empty_seq']}")
    report.append(f"Errors: {stats['error']}")
    
    if results and 'motif_hits' in results[0]:
        hits = [r['motif_hits'] for r in results]
        scores = [r['motif_score'] for r in results]
        report.append(f"\nMotif hits:")
        report.append(f"  Mean: {np.mean(hits):.2f}")
        report.append(f"  Median: {np.median(hits):.0f}")
        report.append(f"  Max: {max(hits)}")
        report.append(f"\nMotif scores:")
        report.append(f"  Mean: {np.mean(scores):.2f}")
        report.append(f"  Min: {min(scores):.2f}")
        report.append(f"  Max: {max(scores):.2f}")
    
    if results and 'prob' in results[0]:
        probs = [r['prob'] for r in results]
        report.append(f"\nProbability scores (no-motif intervals already excluded):")
        report.append(f"  Mean: {np.mean(probs):.2f}")
        report.append(f"  Count: {len(probs)}")

        orientations = [r['orientation'] for r in results]
        from collections import Counter
        counts = Counter(orientations)
        report.append(f"\nOrientation call distribution:")
        for orient, count in counts.most_common():
            pct = 100 * count / len(results)
            report.append(f"  {orient}: {count} ({pct:.1f}%)")
    
    if results and 'motif_orientation' in results[0]:
        orientations = [r['motif_orientation'] for r in results]
        from collections import Counter
        counts = Counter(orientations)
        report.append(f"\nOrientation distribution:")
        for orient, count in counts.most_common():
            pct = 100 * count / len(results)
            report.append(f"  {orient}: {count} ({pct:.1f}%)")
    
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

# ---- Main ----
def main():
    parser = argparse.ArgumentParser(
        description="BED motif finder - scan single intervals for motif presence"
    )
    parser.add_argument("-i", "--input", required=True, help="Input BED file")
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA file")
    parser.add_argument("-m", "--motif", required=True, help="Motif PFM file")
    parser.add_argument("-v", "--vcf", help="Optional VCF file with SNPs")
    parser.add_argument("-o", "--output", default="output_motifs.bed", 
                       help="Output BED file")
    parser.add_argument("--json", help="Output JSON file")
    parser.add_argument("--report", help="Output report file")
    parser.add_argument("--threshold", type=float, default=7.0, 
                       help="PSSM threshold (default=7.0)")
    parser.add_argument("--prob", action="store_true", 
                       help="Probabilistic mode: return orientation probability (0-1) + orientation call (>, <, .)")
    parser.add_argument("--stats", action="store_true", 
                       help="Stats mode: return orientation string (> or <)")
    parser.add_argument("--batch", action="store_true", 
                       help="Batch mode: stream genome (low memory)")
    parser.add_argument("--verbose", action="store_true", 
                       help="Verbose logging")
    
    args = parser.parse_args()
    
    # Setup
    global logger
    logger = setup_logging(args.verbose)
    
    try:
        # Load data
        logger.info("=" * 60)
        intervals = load_bed(args.input)
        genome = load_genome(args.genome, batch_mode=args.batch)
        pssm = load_motif_pfm(args.motif)
        snps = load_snps(args.vcf) if args.vcf else None
        
        # Process
        logger.info("=" * 60)
        start = time.time()
        results, stats = process_bed(
            intervals, genome, pssm, snps,
            threshold=args.threshold,
            probabilistic=args.prob,
            stats_mode=args.stats
        )
        elapsed = time.time() - start
        
        # Save outputs
        logger.info("=" * 60)
        mode = 'probabilistic' if args.prob else ('stats' if args.stats else 'default')
        save_bed(args.output, results, mode=mode)
        
        if args.json:
            save_json(args.json, results)
        
        # Report
        generate_report(results, stats, args.report)
        
        logger.info(f"Completed in {elapsed:.2f}s")
        logger.info(f"Results: {len(results)} intervals with motifs")
        
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=args.verbose)
        sys.exit(1)

if __name__ == "__main__":
    main()