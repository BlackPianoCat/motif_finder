# BEDPE motif finder (no OOP, fully self-contained, with comments)

import argparse
import time
from Bio import SeqIO, motifs
from sortedcontainers import SortedList
import numpy as np
import vcf
from tqdm import tqdm  # for progress bar

# ---- Utility functions ----

def order_key(x):
    # Key for sorting interactions by genomic position
    return (x['chr1'], x['pos1'], x['end1'], x['chr2'], x['pos2'], x['end2'])

def load_interactions(file_path, max_length=0, min_pet=0):
    # Load BEDPE-format interaction file into list of dictionaries
    interactions = SortedList(key=order_key)
    with open(file_path) as f:
        for line in f:
            vals = line.strip().split('\t')
            if max_length > 0 and int(vals[5]) - int(vals[1]) > max_length:
                continue
            chr1 = vals[0] if vals[0].startswith("chr") else "chr" + vals[0]
            chr2 = vals[3] if vals[3].startswith("chr") else "chr" + vals[3]
            pet = int(vals[6]) if vals[6] != "." else 0
            if pet < min_pet:
                continue
            interaction = {
                'chr1': chr1,
                'pos1': int(vals[1]),
                'end1': int(vals[2]),
                'chr2': chr2,
                'pos2': int(vals[4]),
                'end2': int(vals[5]),
                'pet': pet
            }
            interactions.add(interaction)
    print(f"Loaded {len(interactions)} interactions from {file_path}")
    return interactions

def save_interactions(filename, interactions, add_prob=False, add_type=False):
    # Save interactions back to BEDPE with optional prob/type fields
    with open(filename, "w") as out:
        for x in interactions:
            line = f"{x['chr1']}\t{x['pos1']}\t{x['end1']}\t{x['chr2']}\t{x['pos2']}\t{x['end2']}\t{x['pet']}"
            if add_prob:
                line += f"\t{x.get('prob1', '.')}\t{x.get('prob2', '.')}"
            if add_type:
                line += f"\t{x.get('type', '.')}"
            out.write(line + "\n")
    print(f"Saved {len(interactions)} interactions to {filename}")

def load_genome(fasta_path):
    # Load genome FASTA into dictionary of chrom:sequence
    genome = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        if "_" in record.id or "*" in record.id or "EBV" in record.id:
            continue
        genome[record.id] = record.seq
    print(f"Loaded genome with {len(genome)} contigs from {fasta_path}")
    return genome

def load_snps(vcf_path):
    # Load SNPs from VCF into a list of tuples (chr, pos, [allele1, allele2])
    snps = SortedList()
    reader = vcf.Reader(filename=vcf_path)
    count = 0
    for record in reader:
        if "chr" not in record.CHROM or any(x in record.CHROM for x in ["_", "*", "EBV"]) or record.is_indel:
            continue
        chr = record.CHROM
        pos = record.POS
        try:
            # Handle bi-allelic or high/low diversity SNPs
            alleles = [record.ALT[0].sequence, record.ALT[1].sequence] if len(record.ALT) > 1 else \
                      [record.REF, record.ALT[0].sequence] if record.nucl_diversity > 0.99 else \
                      [record.ALT[0].sequence] * 2
            snps.add((chr, pos, alleles))
            count += 1
        except:
            continue
    print(f"Loaded {count} SNPs from {vcf_path}")
    return snps

def update_sequence(seq, snps, offset, allele_id):
    # Modify reference sequence with alleles from SNPs
    for chr, snp_pos, alleles in snps:
        rel_pos = snp_pos - offset
        if 0 <= rel_pos < len(seq):
            seq = seq[:rel_pos] + alleles[allele_id] + seq[rel_pos+1:]
    return seq

# ---- Motif logic ----

def scan_motif(pssm, seq1, seq2, probabilistic=False, get_stats=False):
    # Scan two sequences with PSSM and decide motif presence or orientation
    if probabilistic:
        def side_strength(results):
            left = sum(2**score for orient, score in results if orient < 0)
            right = sum(2**score for orient, score in results if orient >= 0)
            total = left + right
            return -1 if total == 0 else left / total

        r1 = list(pssm.search(seq1, threshold=2.0))
        r2 = list(pssm.search(seq2, threshold=2.0))
        return (side_strength(r1), side_strength(r2))

    r1 = list(pssm.search(seq1, threshold=7.0))
    r2 = list(pssm.search(seq2, threshold=7.0))

    if get_stats:
        orient1 = max(r1, key=lambda x: x[1], default=(0, -1))[0]
        orient2 = max(r2, key=lambda x: x[1], default=(0, -1))[0]
        if not r1 and not r2:
            return "None"
        if not r2:
            return ">.." if orient1 >= 0 else "<.."
        if not r1:
            return "..>" if orient2 >= 0 else "..<"
        if orient1 >= 0 and orient2 < 0:
            return ">..<"
        if orient1 * orient2 >= 0:
            return ">..>" if orient1 >= 0 else "<..<"
        return "<..>"

    for x in r1:
        for y in r2:
            if x[0] >= 0 and y[0] < 0:
                return True
            if x[1] >= 9.0 and y[1] >= 9.0 and x[0]*y[0] >= 0:
                return True
            if x[1] >= 12.0 and y[1] >= 12.0 and x[0] < 0 and y[0] >= 0:
                return True
    return False

def process_interactions(interactions, genome, pssm, snps=None, prob=False, stats=False, ext=0):
    # For each interaction, extract sequences, optionally apply SNPs, and scan motif
    results = SortedList(key=order_key)
    for x in tqdm(interactions, desc="Scanning interactions", ncols=80, colour="cyan"):
        seq1 = genome[x['chr1']][x['pos1']-ext:x['end1']+ext]
        seq2 = genome[x['chr2']][x['pos2']-ext:x['end2']+ext]

        snps1 = [s for s in snps if s[0] == x['chr1'] and x['pos1'] <= s[1] <= x['end1']] if snps else None
        snps2 = [s for s in snps if s[0] == x['chr2'] and x['pos2'] <= s[1] <= x['end2']] if snps else None

        if snps:
            seq1_a = update_sequence(seq1._data, snps1, x['pos1']+1, 0)
            seq2_a = update_sequence(seq2._data, snps2, x['pos2']+1, 0)
            d = scan_motif(pssm, seq1_a, seq2_a, prob, stats)
        else:
            d = scan_motif(pssm, seq1, seq2, prob, stats)

        if prob:
            x['prob1'], x['prob2'] = d
            results.add(x)
        elif stats:
            x['type'] = d
            results.add(x)
        elif d:
            results.add(x)

    return results

# ---- Main entrypoint ----

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input BEDPE file")
    parser.add_argument("-g", "--genome", default="reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa", help="Genome FASTA")
    parser.add_argument("-m", "--motif", default="MA0139.1.pfm", help="Motif PFM file")
    parser.add_argument("-v", "--vcf", help="Optional VCF file with SNPs")
    parser.add_argument("-o", "--output", default="output_with_motifs.bedpe", help="Output BEDPE with motifs")
    parser.add_argument("--prob", action="store_true", help="Enable probabilistic scoring")
    parser.add_argument("--stats", action="store_true", help="Output orientation/stats instead of binary")
    parser.add_argument("--min-pet", type=int, default=0, help="Minimum PET count to keep a loop (default=0)")
    args = parser.parse_args()

    # Load data
    interactions = load_interactions(args.input, min_pet=args.min_pet)
    genome = load_genome(args.genome)
    snps = load_snps(args.vcf) if args.vcf else None

    # Load motif PFM and convert to PSSM
    print(f"Loading motif from {args.motif}")
    with open(args.motif) as f:
        pfm = motifs.read(f, "pfm")
        pssm = pfm.pssm

    # Process all interactions
    print("Scanning motifs...")
    start = time.time()
    result = process_interactions(interactions, genome, pssm, snps, prob=args.prob, stats=args.stats)
    save_interactions(args.output, result, add_prob=args.prob, add_type=args.stats)
    print(f"Done. {len(result)} interactions processed in {time.time() - start:.2f}s.")

if __name__ == "__main__":
    main()
