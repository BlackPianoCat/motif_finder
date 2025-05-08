# Motif Scanning Scripts for Chromatin Interactions and Peaks

This repository contains two standalone Python scripts for scanning DNA motifs in:

1. **3D chromatin interaction regions** (e.g., from BEDPE files)
2. **ChIP-seq or ATAC-seq peak regions** (BED format)

Both tools use Biopython's `motifs` and `pssm` utilities to detect motif presence, orientation, or strand bias.

---

## 🔬 1. `motif_finder_bedpe.py`

Scan motifs in chromatin interaction anchor regions using reference genome and optionally SNPs from a VCF file.

### ✅ Features:

* Input: BEDPE file (e.g. from ChIA-PET, HiChIP)
* Motif scanning using Position-Specific Scoring Matrix (PSSM)
* Supports probabilistic orientation scoring
* Optionally integrates SNP alleles from a VCF file
* Outputs filtered or annotated BEDPE

### 📥 Arguments:

```bash
-i / --input      # Input BEDPE file (required)
-g / --genome     # Reference genome FASTA (default: GRCh38 full set)
-m / --motif      # Motif PFM file (default: MA0139.1.pfm, i.e., CTCF)
-v / --vcf        # Optional VCF file with SNPs
-o / --output     # Output BEDPE file (default: output_with_motifs.bedpe)
--prob            # Flag for probabilistic scoring
--stats           # Flag to output motif orientation category
--min-pet         # To filter out loops with value less than this strength
```

To download the fasta files: https://knowledge.illumina.com/software/on-premises-software/software-on-premises-software-troubleshooting-list/000007409.

### 🚀 Example:

```bash
python motif_finder_bedpe.py \
  -i loops.bedpe \
  -g GRCh38.fa \
  -m MA0139.1.pfm \
  -v sample.vcf \
  -o loops_filtered.bedpe \
  --prob
```

---

## 🧬 2. `chipseq_motif_finder.py`

Scan motifs over peak regions from ChIP-seq or ATAC-seq.

### ✅ Features:

* Input: BED/narrowPeak file
* Extracts sequences from reference genome
* Scans with motif PSSM (threshold-based)
* Outputs hits in BED format

### 📥 Arguments:

```bash
-i / --input      # Input BED file with peaks
-g / --genome     # Genome FASTA file
-m / --motif      # Motif PFM file (default: MA0139.1.pfm)
-o / --output     # Output BED file for motif hits
--prob            # Enable probabilistic mode (outputs 1bp hits with scores)
-t / --threshold  # PSSM match threshold (default: 7.0)
```

### 🚀 Example:

```bash
python chipseq_motif_finder.py \
  -i peaks.bed \
  -g GRCh38.fa \
  -m MA0139.1.pfm \
  -o ctcf_hits.bed \
  --threshold 6.5
```

---

## 🔧 Requirements:

```bash
pip install biopython sortedcontainers vcf tqdm
```

---

## 📁 Input Format Notes:

* Genome FASTA should use `chr` names matching BED/BEDPE
* Motif PFM files should be in JASPAR/TRANSFAC-compatible format
* VCF file should contain bi-allelic SNPs (INDELs are ignored)

---

## 📜 License

MIT License — free to use and adapt. Attribution appreciated!
Feel free to request enhancements or extra features.
