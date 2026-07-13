# CTCF Motif Finder Suite

Three related scripts for scanning genomic regions for transcription-factor motifs (e.g. CTCF) using a position-specific scoring matrix (PSSM). They share the same core scanning engine and three output modes, but differ in what kind of input region they expect.

| Script | Input format | Region shape |
|---|---|---|
| `motif_finder_bedpe.py` | BEDPE (ChIA-PET / HiChIP loop calls) | **Pair** of anchors (chr1/pos1/end1 <-> chr2/pos2/end2) |
| `motif_finder_peaks.py` | narrowPeak (MACS2 / HiChIP peaks) | **Single** peak with a summit |
| `motif_finder_chipseq.py` | BED (any generic interval) | **Single** plain interval |

All three:
- Load a genome FASTA (in-memory or `--batch` streaming mode for low memory use)
- Load one (or, for the BEDPE version, multiple) motif PFM file(s) and build a PSSM with Biopython
- Optionally inject SNPs from a VCF into the sequence before scanning
- Support the same three scanning modes via CLI flags
- Emit a text output, optional JSON, and a summary report

---

## Core algorithmic logic

Each region's sequence is scanned with `pssm.search(seq, threshold=...)`, which returns `(position, score)` hits. A **positive position** means the motif matched the forward strand; a **negative position** means it matched the reverse complement. This single primitive powers all three modes.

### The underlying math (PSSM search)

The PFM file gives a base count at each of the motif's *L* positions. Biopython converts this into a **position-specific scoring matrix (PSSM)** of log-odds values:

```
PSSM[b, i] = log2( P(base=b at position i) / P(base=b background) )
```

where `P(base=b at position i)` comes from the normalized PFM counts (with a pseudocount) and the background is typically a uniform 25% per base. `pssm.search` then slides this *L*-bp window across the sequence, forward strand and reverse complement both, and at each position `j` sums the matrix entries for the bases actually observed:

```
score(j) = sum_i PSSM[seq[j+i], i]   for i = 0 .. L-1
```

A `score(j)` above the chosen `--threshold` is reported as a hit at position `j` (forward) or `-j` (reverse complement), so higher scores mean a better match to the motif's consensus and the threshold is effectively a minimum log-odds bar a window must clear to count as a real binding site. Everything downstream (hit counts, best score, `2**score`-weighted orientation bias, orientation strings) is built directly from these `(position, score)` pairs.

1. **Binary mode (default)** - keep the region only if at least one hit clears the threshold. Extra metrics (hit count, best score, orientation, distance to an anchor/summit) are computed for the surviving regions. The BEDPE version additionally checks for a **canonical inverted pair** (forward hit on the left anchor, reverse hit on the right anchor - the classic CTCF loop-anchor convergence signal), which can be enforced with `--canonical-only`.
2. **Probabilistic mode (`--prob`)** - instead of a hard yes/no, compute a continuous **orientation-bias score** between 0 and 1: `sum(2^score for reverse hits) / sum(2^score for all hits)`. Higher-scoring hits dominate the weighting. Returns `-1` if no hits are found at all.
   - BEDPE: bias is computed per anchor and reported as `prob1`/`prob2`.
   - narrowPeak/BED: bias is computed once per single sequence, expressed as forward-vs-reverse instead of left-vs-right.
3. **Stats mode (`--stats`)** - collapse hits into a compact orientation string:
   - BEDPE: `">..<"` (canonical), `"<..>"` (inverted), `">..>"`/`"<..<"` (parallel), or `"None"`.
   - narrowPeak: `">"`, `"<"`, `"><"` (both orientations present), or `"None"`.
   - BED: `">"`, `"<"`, or `"None"` (based on the single best-scoring hit).

The three modes are mutually exclusive at runtime - pass at most one of `--prob` / `--stats`; omitting both runs binary mode.

---

## 1. `motif_finder_bedpe.py`

Scans **both anchors** of a chromatin loop call and evaluates whether the motif pair supports a CTCF-style convergent orientation.

**Extra features:** multi-motif support (`-m motif1.pfm motif2.pfm ...`), distance-to-anchor filtering, quality score (average motif score + canonical-orientation bonus).

```bash
python3 motif_finder_bedpe.py \
  -i loops.bedpe -g genome.fa -m CTCF.pfm \
  -o out.bedpe --canonical-only --max-distance 500 \
  --include-scores --json out.json --report report.txt
```

Key flags: `--canonical-only`, `--max-distance`, `--min-pet`, `--include-scores`, `--prob`, `--stats`, `--batch`.

---

## 2. `motif_finder_peaks.py`

Scans a **single narrowPeak region** and reports orientation/quality relative to the peak **summit** (column 10), since there's no second anchor to pair against.

**Extra features:** forward/reverse hit counts, distance from best hit to summit, quality bonus if the hit sits within 20 bp of the summit, `--min-score` / `--min-qvalue` filters on the narrowPeak columns.

```bash
python3 motif_finder_peaks.py \
  -i peaks.narrowPeak -g genome.fa -m CTCF.pfm \
  -o out.narrowPeak --max-distance 300 --min-qvalue 2 \
  --include-scores --json out.json --report report.txt
```

Key flags: `--min-score`, `--min-qvalue`, `--max-distance`, `--include-scores`, `--prob`, `--stats`, `--batch`.

---

## 3. `motif_finder_chipseq.py`

The simplest variant - scans a **plain BED interval** with no anchor pairing and no summit; just "does this interval contain the motif, and with what orientation/score."

**Extra features:** SNP-aware sequence editing via `--vcf`, minimal output (hit count + best score in default mode).

```bash
python3 motif_finder_chipseq.py \
  -i regions.bed -g genome.fa -m CTCF.pfm \
  -v snps.vcf -o out.bed --json out.json --report report.txt
```

Key flags: `--vcf`, `--threshold`, `--prob`, `--stats`, `--batch`.

---

## Common requirements

```bash
pip install biopython sortedcontainers numpy tqdm PyVCF
```

- Motif files must be in **PFM (position frequency matrix)** format readable by `Bio.motifs.read(..., "pfm")`.
- Chromosome names are auto-normalized to the `chr1` style (a bare `1` becomes `chr1`).
- `--batch` trades speed for memory by re-parsing the FASTA per region instead of loading the whole genome into RAM - use it for large/many-contig genomes.