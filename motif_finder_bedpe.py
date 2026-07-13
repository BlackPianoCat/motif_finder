#!/usr/bin/env python3
# ============================================================================
# ENHANCED BEDPE motif finder with robust error handling & cool features
# ============================================================================

import argparse
import time
import json
import os
import sys
import logging
from pathlib import Path
from Bio import SeqIO, motifs
from sortedcontainers import SortedList
import numpy as np
import vcf
from tqdm import tqdm
from collections import defaultdict, Counter
from dataclasses import dataclass, asdict

# ---- Setup logging ----
def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
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

@dataclass
class InteractionResult:
    """Enhanced interaction with motif results"""
    chr1: str
    pos1: int
    end1: int
    chr2: str
    pos2: int
    end2: int
    pet: int
    motif1_hits: list = None  # List of MotifHit
    motif2_hits: list = None  # List of MotifHit
    canonical: bool = False  # True if inverted pair (CTCF-like)
    prob1: float = None
    prob2: float = None
    dist_to_motif1: int = None
    dist_to_motif2: int = None
    quality_score: float = None

# ---- Utility functions ----
def order_key(x):
    """Key for sorting interactions by genomic position"""
    if isinstance(x, dict):
        return (x['chr1'], x['pos1'], x['end1'], x['chr2'], x['pos2'], x['end2'])
    else:  # InteractionResult
        return (x.chr1, x.pos1, x.end1, x.chr2, x.pos2, x.end2)

def file_exists(file_path, context=""):
    """Validate file existence"""
    if not Path(file_path).exists():
        raise FileNotFoundError(f"{context}: {file_path} not found")
    return True

def validate_bedpe_format(line, line_num):
    """Validate BEDPE line format"""
    vals = line.strip().split('\t')
    if len(vals) < 7:
        raise ValueError(f"Line {line_num}: Expected ≥7 BEDPE fields, got {len(vals)}")
    try:
        int(vals[1]), int(vals[2]), int(vals[4]), int(vals[5])
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
def load_interactions(file_path, max_length=0, min_pet=0):
    """Load BEDPE-format interaction file with validation"""
    logger.info(f"Loading interactions from {file_path}")
    file_exists(file_path, "BEDPE file")
    
    interactions = SortedList(key=order_key)
    skipped = 0
    total_lines = 0
    blank_or_comment = 0
    
    try:
        with open(file_path) as f:
            for line_num, line in enumerate(f, 1):
                total_lines += 1
                if line.startswith("#") or not line.strip():
                    blank_or_comment += 1
                    continue
                try:
                    vals = validate_bedpe_format(line, line_num)
                    
                    chr1 = normalize_chrom(vals[0])
                    chr2 = normalize_chrom(vals[3])
                    pos1, end1 = int(vals[1]), int(vals[2])
                    pos2, end2 = int(vals[4]), int(vals[5])
                    pet = int(vals[6]) if vals[6] != "." else 0
                    
                    # Validation checks
                    if max_length > 0 and max(end1 - pos1, end2 - pos2) > max_length:
                        skipped += 1
                        continue
                    if pet < min_pet:
                        skipped += 1
                        continue
                    
                    interaction = {
                        'chr1': chr1, 'pos1': pos1, 'end1': end1,
                        'chr2': chr2, 'pos2': pos2, 'end2': end2,
                        'pet': pet
                    }
                    interactions.add(interaction)
                except ValueError as e:
                    logger.warning(f"Skipping line {line_num}: {e}")
                    skipped += 1
    except IOError as e:
        raise IOError(f"Cannot read BEDPE file: {e}")
    
    logger.info(
        f"Loaded {len(interactions)} interactions (skipped {skipped} malformed, "
        f"{blank_or_comment} blank/comment, {total_lines} lines total in file)"
    )
    return interactions

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
            # Filter out contaminants
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

def load_snps(vcf_path):
    """Load SNPs from VCF with quality filtering"""
    file_exists(vcf_path, "VCF file")
    logger.info(f"Loading SNPs from {vcf_path}")
    
    snps = SortedList()
    count = 0
    
    try:
        reader = vcf.Reader(filename=vcf_path)
        for record in reader:
            # Filter criteria
            if not record.CHROM or record.is_indel or record.is_sv:
                continue
            if any(x in record.CHROM for x in ["_", "*", "EBV", "decoy"]):
                continue
            
            chr_name = normalize_chrom(record.CHROM)
            pos = record.POS
            
            try:
                # Extract alleles
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

def update_sequence(seq, snps, offset):
    """Apply SNP alleles to sequence"""
    seq_list = list(str(seq))
    for chr_name, snp_pos, alleles in snps:
        rel_pos = snp_pos - offset
        if 0 <= rel_pos < len(seq_list):
            seq_list[rel_pos] = alleles[0]  # Use first allele
    return ''.join(seq_list)

# ---- Motif scanning (restored original logic with enhancements) ----
def scan_sequence_raw(pssm, seq, threshold=7.0):
    """Scan sequence and return raw pssm.search results (position, score)"""
    try:
        return list(pssm.search(seq, threshold=threshold))
    except Exception as e:
        logger.debug(f"Error scanning sequence: {e}")
        return []

def scan_sequence(pssm, seq, threshold=7.0):
    """Scan sequence and return all hits above threshold"""
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

def side_strength(results):
    """
    RESTORED FROM ORIGINAL: Calculate probabilistic orientation
    Returns probability that motif is on left (negative orientation)
    Returns -1 if no motifs found
    """
    if not results:
        return -1.0
    
    left = sum(2**score for orient, score in results if orient < 0)
    right = sum(2**score for orient, score in results if orient >= 0)
    total = left + right
    
    if total == 0:
        return -1.0
    return left / total

def scan_motif_probabilistic(pssm, seq1, seq2, threshold=2.0):
    """
    RESTORED FROM ORIGINAL: Probabilistic scoring mode
    Returns tuple of probabilities (prob1, prob2)
    Each probability: 0-1 (left-biased to right-biased) or -1 (no motifs)
    """
    r1 = scan_sequence_raw(pssm, seq1, threshold=threshold)
    r2 = scan_sequence_raw(pssm, seq2, threshold=threshold)
    
    prob1 = side_strength(r1)
    prob2 = side_strength(r2)
    
    return (prob1, prob2)

def scan_motif_stats(pssm, seq1, seq2, threshold=7.0):
    """
    RESTORED FROM ORIGINAL: Stats/orientation string mode
    Returns string describing orientation pattern:
    - ">...<" : forward on left, reverse on right (canonical CTCF)
    - ">.." : only left has forward
    - "..<" : only right has reverse
    - "<..>" : reverse on left, forward on right (inverted)
    - ">..>" : both forward
    - "<..<" : both reverse
    - "None" : no motifs found
    """
    r1 = scan_sequence_raw(pssm, seq1, threshold=threshold)
    r2 = scan_sequence_raw(pssm, seq2, threshold=threshold)
    
    # Get orientation of best hit (if any)
    orient1 = max(r1, key=lambda x: x[1], default=(0, -1))[0]
    orient2 = max(r2, key=lambda x: x[1], default=(0, -1))[0]
    
    # No motifs found
    if not r1 and not r2:
        return "None"
    
    # Only one side has motif
    if not r2:
        return ">.." if orient1 >= 0 else "< .."
    if not r1:
        return "..>" if orient2 >= 0 else "..< "
    
    # Both sides have motifs - determine configuration
    if orient1 >= 0 and orient2 < 0:
        return ">..<"  # Canonical
    if orient1 * orient2 >= 0:
        return ">..>" if orient1 >= 0 else "<..<"  # Parallel
    else:
        return "<..>"  # Inverted

def scan_motif_binary(pssm, seq1, seq2, threshold=7.0):
    """
    RESTORED FROM ORIGINAL: Binary mode with original heuristics
    Returns True if motif pair matches biological criteria, False otherwise
    
    Criteria:
    1. One forward, one reverse with high scores (canonical CTCF)
    2. Both high-scoring and same orientation (cooperative binding)
    3. Very high scores with opposite orientations
    """
    r1 = scan_sequence_raw(pssm, seq1, threshold=threshold)
    r2 = scan_sequence_raw(pssm, seq2, threshold=threshold)
    
    # Original heuristic logic
    for x in r1:
        for y in r2:
            # Canonical: forward on left, reverse on right
            if x[0] >= 0 and y[0] < 0:
                return True
            # Both high-scoring and same orientation
            if x[1] >= 9.0 and y[1] >= 9.0 and x[0]*y[0] >= 0:
                return True
            # Very high scores with specific orientation
            if x[1] >= 12.0 and y[1] >= 12.0 and x[0] < 0 and y[0] >= 0:
                return True
    return False

def detect_canonical_pair(hits1, hits2):
    """
    Detect canonical CTCF orientation (inverted pair: > ... <)
    CTCF has a strong preferential binding orientation for chromatin looping
    """
    if not hits1 or not hits2:
        return False, None
    
    # Best hits
    best1 = max(hits1, key=lambda h: h.score)
    best2 = max(hits2, key=lambda h: h.score)
    
    # Canonical: forward on left, reverse on right
    is_canonical = best1.orientation > 0 and best2.orientation < 0
    orientation_pair = (best1.orientation, best2.orientation)
    
    return is_canonical, orientation_pair

def orientation_from_hits(hits1, hits2):
    """
    Build the same compact arrow string as scan_motif_stats() (">..<", ">..",
    "..<", "<..>", ">..>", "<..<", or "None"), but from MotifHit lists that
    were already scanned in default/binary mode - so default mode can show
    the CTCF orientation without having to re-scan in --stats mode.
    """
    best1 = max(hits1, key=lambda h: h.score) if hits1 else None
    best2 = max(hits2, key=lambda h: h.score) if hits2 else None

    if best1 is None and best2 is None:
        return "None"
    if best1 is None:
        return "..>" if best2.orientation >= 0 else "..<"
    if best2 is None:
        return ">.." if best1.orientation >= 0 else "<.."

    if best1.orientation >= 0 and best2.orientation < 0:
        return ">..<"  # Canonical
    if best1.orientation * best2.orientation >= 0:
        return ">..>" if best1.orientation >= 0 else "<..<"  # Parallel
    return "<..>"  # Inverted

def calculate_motif_distance(hits, anchor_start, anchor_end):
    """Calculate minimum distance from motif hit to anchor region"""
    if not hits:
        return None
    
    best_hit = max(hits, key=lambda h: h.score)
    
    # Distance to nearest edge of anchor
    if best_hit.position < anchor_start:
        return anchor_start - best_hit.position
    elif best_hit.position > anchor_end:
        return best_hit.position - anchor_end
    else:
        return 0  # Motif is within anchor

def calculate_quality_score(hits1, hits2, canonical=False, distance_weight=0.1):
    """
    Calculate overall quality score combining:
    - Motif scores
    - Canonical orientation
    - Distance to anchor (bonus for close proximity)
    """
    if not hits1 and not hits2:
        return 0.0
    
    score = 0.0
    
    if hits1:
        avg_score1 = np.mean([h.score for h in hits1])
        score += avg_score1 / 2
    
    if hits2:
        avg_score2 = np.mean([h.score for h in hits2])
        score += avg_score2 / 2
    
    if canonical:
        score += 2.0  # Bonus for canonical orientation
    
    return min(score, 20.0)  # Cap at 20

# ---- Main processing ----
def process_interactions(interactions, genome, pssms, snps=None, threshold=7.0, 
                        distance_threshold=None, canonical_only=False,
                        probabilistic=False, stats_mode=False):
    """
    Process all interactions and scan for motifs
    
    Modes:
    - Default: Return interactions with motifs (binary True/False)
    - probabilistic=True: Return probability scores (-1 to 1 for each anchor)
    - stats_mode=True: Return orientation strings (">...<", "None", etc.)
    """
    logger.info(f"Scanning {len(interactions)} interactions with {len(pssms)} motif(s)")
    if probabilistic:
        logger.info("Using probabilistic scoring mode")
    if stats_mode:
        logger.info("Using stats/orientation mode")
    
    results = SortedList(key=order_key)
    stats = defaultdict(int)
    
    for interaction in tqdm(interactions, desc="Scanning", ncols=90, colour="cyan"):
        try:
            chr1, pos1, end1 = interaction['chr1'], interaction['pos1'], interaction['end1']
            chr2, pos2, end2 = interaction['chr2'], interaction['pos2'], interaction['end2']
            
            # Extract sequences
            if isinstance(genome, str):  # Batch mode
                seq1 = get_sequence_batch(genome, chr1, pos1, end1)
                seq2 = get_sequence_batch(genome, chr2, pos2, end2)
            else:
                seq1 = genome.get(chr1, None)
                seq2 = genome.get(chr2, None)
            
            if seq1 is None or seq2 is None:
                logger.debug(f"Chromosome not found: {chr1} or {chr2}")
                stats['missing_chrom'] += 1
                continue
            
            # Get subsequences
            seq1 = seq1[pos1:end1]
            seq2 = seq2[pos2:end2]
            
            if len(seq1) == 0 or len(seq2) == 0:
                stats['empty_seq'] += 1
                continue
            
            # Apply SNPs if available
            if snps:
                snps1 = [s for s in snps if s[0] == chr1 and pos1 <= s[1] <= end1]
                snps2 = [s for s in snps if s[0] == chr2 and pos2 <= s[1] <= end2]
                if snps1:
                    seq1 = update_sequence(seq1, snps1, pos1)
                if snps2:
                    seq2 = update_sequence(seq2, snps2, pos2)
            
            seq1_str = str(seq1)
            seq2_str = str(seq2)
            
            # Process based on mode
            if probabilistic:
                # RESTORED: Probabilistic mode - returns probability scores
                for motif_id, pssm in pssms.items():
                    prob1, prob2 = scan_motif_probabilistic(pssm, seq1_str, seq2_str, threshold=2.0)
                    
                    # Store probabilities (-1 if no motifs, 0-1 otherwise)
                    result = {
                        'chr1': chr1, 'pos1': pos1, 'end1': end1,
                        'chr2': chr2, 'pos2': pos2, 'end2': end2,
                        'pet': interaction['pet'],
                        'prob1': prob1,  # -1 if no motif, else 0-1
                        'prob2': prob2,  # -1 if no motif, else 0-1
                        'motif_id': motif_id
                    }
                    results.add(result)
                    stats['probabilistic'] += 1
                    
            elif stats_mode:
                # RESTORED: Stats mode - returns orientation strings
                for motif_id, pssm in pssms.items():
                    orientation_type = scan_motif_stats(pssm, seq1_str, seq2_str, threshold)
                    
                    # Store orientation pattern
                    result = {
                        'chr1': chr1, 'pos1': pos1, 'end1': end1,
                        'chr2': chr2, 'pos2': pos2, 'end2': end2,
                        'pet': interaction['pet'],
                        'type': orientation_type,  # ">..<", "None", etc.
                        'motif_id': motif_id
                    }
                    results.add(result)
                    stats['stats'] += 1
                    
            else:
                # Default: Binary mode with enhanced metrics
                # Scan all motifs
                all_hits1 = []
                all_hits2 = []
                passed_filter = False
                
                for motif_id, pssm in pssms.items():
                    # Use binary mode logic from original
                    has_motif = scan_motif_binary(pssm, seq1_str, seq2_str, threshold)
                    
                    if has_motif:
                        passed_filter = True
                        # Also get detailed hits for metrics
                        hits1 = scan_sequence(pssm, seq1_str, threshold)
                        hits2 = scan_sequence(pssm, seq2_str, threshold)
                        
                        for h in hits1:
                            h.motif_id = motif_id
                        for h in hits2:
                            h.motif_id = motif_id
                        
                        all_hits1.extend(hits1)
                        all_hits2.extend(hits2)
                
                if not passed_filter:
                    stats['no_motif'] += 1
                    continue
                
                # Check canonical
                is_canonical, _ = detect_canonical_pair(all_hits1, all_hits2)
                if canonical_only and not is_canonical:
                    stats['not_canonical'] += 1
                    continue

                # Arrow orientation string (">..<" etc.), always computed so
                # default-mode output can show CTCF orientation without --stats
                orientation_type = orientation_from_hits(all_hits1, all_hits2)

                # Calculate metrics
                dist1 = calculate_motif_distance(all_hits1, pos1, end1)
                dist2 = calculate_motif_distance(all_hits2, pos2, end2)
                quality = calculate_quality_score(all_hits1, all_hits2, is_canonical)
                
                # Distance filtering
                if distance_threshold and ((dist1 and dist1 > distance_threshold) or 
                                         (dist2 and dist2 > distance_threshold)):
                    stats['far_from_anchor'] += 1
                    continue
                
                # Create result
                result = {
                    'chr1': chr1, 'pos1': pos1, 'end1': end1,
                    'chr2': chr2, 'pos2': pos2, 'end2': end2,
                    'pet': interaction['pet'],
                    'motif_hits1': len(all_hits1),
                    'motif_hits2': len(all_hits2),
                    'canonical': is_canonical,
                    'orientation': orientation_type,  # ">..<", "<..>", "None", etc.
                    'dist_motif1': dist1,
                    'dist_motif2': dist2,
                    'quality_score': quality
                }
                
                results.add(result)
                stats['passed'] += 1
            
        except Exception as e:
            logger.debug(f"Error processing interaction {interaction}: {e}")
            stats['error'] += 1
    
    return results, stats

# ---- Output functions ----
def save_bedpe(filename, results, include_scores=False):
    """
    Save results to BEDPE format
    Handles three output formats:
    - Default: BEDPE + orientation arrow (">..<", "<..>", "None", etc.), always
      shown; --include-scores additionally appends motif count/quality/canonical/distance
    - Probabilistic (--prob): BEDPE + prob1 + prob2 (each -1 if no motifs, else 0-1)
    - Stats (--stats): BEDPE + orientation type string
    """
    logger.info(f"Writing {len(results)} results to {filename}")
    
    try:
        with open(filename, "w") as f:
            for r in results:
                line = f"{r['chr1']}\t{r['pos1']}\t{r['end1']}\t{r['chr2']}\t{r['pos2']}\t{r['end2']}\t{r['pet']}"
                
                # Probabilistic mode: add probability scores
                if 'prob1' in r and 'prob2' in r:
                    prob1_str = f"{r['prob1']:.2f}" if r['prob1'] != -1 else "-1"
                    prob2_str = f"{r['prob2']:.2f}" if r['prob2'] != -1 else "-1"
                    line += f"\t{prob1_str}\t{prob2_str}"
                    if 'motif_id' in r:
                        line += f"\t{r['motif_id']}"
                
                # Stats mode: add orientation type
                elif 'type' in r:
                    line += f"\t{r['type']}"
                    if 'motif_id' in r:
                        line += f"\t{r['motif_id']}"
                
                # Default mode: always show the orientation arrow;
                # --include-scores additionally appends the numeric detail columns
                elif 'motif_hits1' in r:
                    line += f"\t{r['orientation']}"
                    if include_scores:
                        line += f"\t{r['motif_hits1']}\t{r['motif_hits2']}\t{r['quality_score']:.2f}"
                        line += f"\t{r['canonical']}\t{r.get('dist_motif1', '.')}\t{r.get('dist_motif2', '.')}"
                
                f.write(line + "\n")
    except IOError as e:
        raise IOError(f"Cannot write output file: {e}")

def save_json(filename, results, full_results=True):
    """Save results as JSON for downstream analysis"""
    logger.info(f"Writing JSON results to {filename}")
    
    try:
        with open(filename, "w") as f:
            json.dump([dict(r) for r in results], f, indent=2, default=str)
    except Exception as e:
        raise IOError(f"Cannot write JSON file: {e}")

def generate_report(results, stats, output_file=None):
    """
    Generate summary statistics report.

    FIX: the report body is now mode-aware. The old version always tried to
    read 'quality_score' / 'motif_hits1' / 'canonical' from every result,
    which only exist in default (binary) mode. In --prob mode results only
    have 'prob1'/'prob2'/'motif_id', and in --stats mode only 'type' - reading
    the wrong keys raised `KeyError: 'quality_score'` (surfaced as
    "Fatal error: 'quality_score'"). We now branch on which keys are actually
    present before computing mode-specific stats.
    """
    report = []
    report.append("\n" + "="*60)
    report.append("MOTIF SCANNING REPORT")
    report.append("="*60)
    report.append(f"\nTotal interactions passed filters: {stats['passed']}")
    report.append(f"Skipped - no motif found: {stats['no_motif']}")
    report.append(f"Skipped - not canonical orientation: {stats['not_canonical']}")
    report.append(f"Skipped - motif too far from anchor: {stats['far_from_anchor']}")
    report.append(f"Skipped - missing chromosome: {stats['missing_chrom']}")
    report.append(f"Skipped - empty sequence: {stats['empty_seq']}")
    report.append(f"Errors encountered: {stats['error']}")

    if results and 'quality_score' in results[0]:
        # ---- Default (binary) mode ----
        canonical_count = sum(1 for r in results if r.get('canonical', False))
        avg_quality = np.mean([r['quality_score'] for r in results])

        report.append(f"\nCanonical pairs: {canonical_count}/{stats['passed']} ({100*canonical_count/stats['passed']:.1f}%)")
        report.append(f"Average quality score: {avg_quality:.2f}")

        # Motif hit distribution
        hits1_dist = [r['motif_hits1'] for r in results]
        hits2_dist = [r['motif_hits2'] for r in results]
        report.append(f"\nMotif hits per anchor:")
        report.append(f"  Anchor 1: mean={np.mean(hits1_dist):.2f}, max={max(hits1_dist)}")
        report.append(f"  Anchor 2: mean={np.mean(hits2_dist):.2f}, max={max(hits2_dist)}")

        # Distance stats
        valid_dists = [d for r in results for d in [r['dist_motif1'], r['dist_motif2']] if d is not None]
        if valid_dists:
            report.append(f"\nDistance to anchor (bp):")
            report.append(f"  Mean: {np.mean(valid_dists):.0f}")
            report.append(f"  Median: {np.median(valid_dists):.0f}")
            report.append(f"  Max: {max(valid_dists):.0f}")

    elif results and 'prob1' in results[0]:
        # ---- Probabilistic mode ----
        probs1 = [r['prob1'] for r in results if r['prob1'] != -1]
        probs2 = [r['prob2'] for r in results if r['prob2'] != -1]
        report.append(f"\nProbabilistic results: {stats['probabilistic']}")
        if probs1:
            report.append(f"  Anchor 1 - mean prob: {np.mean(probs1):.2f}, with motif: {len(probs1)}/{len(results)}")
        if probs2:
            report.append(f"  Anchor 2 - mean prob: {np.mean(probs2):.2f}, with motif: {len(probs2)}/{len(results)}")

    elif results and 'type' in results[0]:
        # ---- Stats mode ----
        types = [r['type'] for r in results]
        counts = Counter(types)
        report.append(f"\nOrientation type distribution:")
        for t, c in counts.most_common():
            pct = 100 * c / len(results)
            report.append(f"  {t}: {c} ({pct:.1f}%)")

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
        description="Enhanced BEDPE motif finder with multi-motif support and robust validation"
    )
    parser.add_argument("-i", "--input", required=True, help="Input BEDPE file")
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA file")
    parser.add_argument("-m", "--motif", required=True, nargs='+', 
                       help="One or more motif PFM files")
    parser.add_argument("-v", "--vcf", help="Optional VCF file with SNPs")
    parser.add_argument("-o", "--output", default="output_motifs.bedpe", 
                       help="Output BEDPE file")
    parser.add_argument("--json", help="Output JSON file with full results")
    parser.add_argument("--report", help="Output statistics report")
    parser.add_argument("--threshold", type=float, default=7.0, 
                       help="PSSM score threshold (default=7.0)")
    parser.add_argument("--canonical-only", action="store_true", 
                       help="Only keep canonical CTCF orientation (> ... <)")
    parser.add_argument("--max-distance", type=int, 
                       help="Max distance from motif to anchor (bp)")
    parser.add_argument("--min-pet", type=int, default=0, 
                       help="Minimum PET count")
    parser.add_argument("--include-scores", action="store_true", 
                       help="Include motif scores in BEDPE output")
    parser.add_argument("--prob", action="store_true", 
                       help="RESTORED: Probabilistic scoring mode (returns -1 to 1 orientation probability)")
    parser.add_argument("--stats", action="store_true", 
                       help="RESTORED: Stats mode (returns orientation strings like >..< or None)")
    parser.add_argument("--batch", action="store_true", 
                       help="Batch mode: stream genome (slower, less memory)")
    parser.add_argument("--verbose", action="store_true", 
                       help="Verbose logging")
    
    args = parser.parse_args()
    
    # Setup logging
    global logger
    logger = setup_logging(args.verbose)
    
    try:
        # Validate inputs
        if isinstance(args.motif, list) and len(args.motif) > 0:
            motif_files = args.motif
        else:
            motif_files = [args.motif]
        
        # Load data
        logger.info("=" * 60)
        interactions = load_interactions(args.input, min_pet=args.min_pet)
        genome = load_genome(args.genome, batch_mode=args.batch)
        pssms = load_multiple_motifs(motif_files)
        snps = load_snps(args.vcf) if args.vcf else None
        
        # Process
        logger.info("=" * 60)
        start = time.time()
        results, stats = process_interactions(
            interactions, genome, pssms, snps,
            threshold=args.threshold,
            distance_threshold=args.max_distance,
            canonical_only=args.canonical_only,
            probabilistic=args.prob,
            stats_mode=args.stats
        )
        elapsed = time.time() - start
        
        # Save outputs
        logger.info("=" * 60)
        save_bedpe(args.output, results, include_scores=args.include_scores)
        if args.json:
            save_json(args.json, results)
        
        # Report
        report_file = args.report if args.report else None
        generate_report(results, stats, report_file)
        
        logger.info(f"Completed in {elapsed:.2f}s")
        logger.info(f"Results saved to {args.output}")
        
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=args.verbose)
        sys.exit(1)

if __name__ == "__main__":
    main()