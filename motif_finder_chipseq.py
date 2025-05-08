# Simplified motif scanner for ChIP/ATAC-seq peak BEDs with progress bars and comments

import argparse
from Bio import SeqIO, motifs
from Bio.Seq import Seq
from tqdm import tqdm
import os
import math

# Load genome FASTA as dictionary {chrom: sequence}
def load_genome(fasta_path):
    genome = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        if not any(x in record.id for x in ["_", "*", "EBV"]):
            genome[record.id] = record.seq
    print(f"Loaded {len(genome)} genome sequences from {fasta_path}")
    return genome

# Load BED/narrowPeak format into list of (chr, start, end) tuples
def parse_bed(bed_file):
    peaks = []
    with open(bed_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            peaks.append((chrom, start, end))
    print(f"Loaded {len(peaks)} peaks from {bed_file}")
    return peaks

# Compute probabilistic score as fraction of signal on the left strand
def compute_probabilistic_score(results):
    left = sum(math.pow(2, score) for orient, score in results if orient < 0)
    right = sum(math.pow(2, score) for orient, score in results if orient >= 0)
    total = left + right
    return -1 if total == 0 else left / total

# Search a single sequence with the motif PSSM, return hits
# If probabilistic=True, return (chr, pos, prob_left)
# Otherwise return BED-like (chr, start, end, score, strand)
def scan_sequence(pssm, sequence, chrom, start, threshold=7.0, probabilistic=False):
    results = list(pssm.search(sequence, threshold=2.0 if probabilistic else threshold))
    if probabilistic:
        prob = compute_probabilistic_score(results)
        return [(chrom, start, start + len(sequence), prob)]
    else:
        hits = []
        for pos, score in results:
            abs_pos = start + pos if pos >= 0 else start + len(sequence) + pos
            strand = "+" if pos >= 0 else "-"
            hits.append((chrom, abs_pos, abs_pos + len(pssm), score, strand))
        return hits

# Main execution logic
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="BED file with ChIP-Seq or ATAC-Seq peaks")
    parser.add_argument("-g", "--genome", required=True, help="Reference genome FASTA")
    parser.add_argument("-m", "--motif", required=True, help="PFM motif file (JASPAR format)")
    parser.add_argument("-o", "--output", required=True, help="Output BED file with motif hits")
    parser.add_argument("--prob", action="store_true", help="Enable probabilistic mode (returns fraction of left-oriented signal)")
    parser.add_argument("-t", "--threshold", type=float, default=7.0, help="PSSM score threshold")
    args = parser.parse_args()

    # Load input data
    genome = load_genome(args.genome)
    peaks = parse_bed(args.input)
    
    # Load motif
    print(f"Loading motif from {args.motif}")
    with open(args.motif) as f:
        pfm = motifs.read(f, "pfm")
        pssm = pfm.pssm

    # Process each peak and collect motif hits
    hits = []
    print("Scanning for motifs across peaks...")
    for chrom, start, end in tqdm(peaks, desc="Scanning peaks", ncols=80, colour="green"):
        if chrom not in genome:
            continue
        sequence = genome[chrom][start:end]
        hits.extend(scan_sequence(pssm, sequence, chrom, start, threshold=args.threshold, probabilistic=args.prob))

    # Write results to output file
    with open(args.output, "w") as out:
        if args.prob:
            for chrom, start, end, prob in hits:
                out.write(f"{chrom}\t{start}\t{end}\t.\t{prob:.3f}\n")
        else:
            for chrom, start, end, score, strand in hits:
                out.write(f"{chrom}\t{start}\t{end}\t.\t{score:.3f}\t{strand}\n")

    print(f"Wrote {len(hits)} motif hits to {args.output}")

if __name__ == "__main__":
    main()
