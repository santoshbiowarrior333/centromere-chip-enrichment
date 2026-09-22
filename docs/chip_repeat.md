# chip_repeat: ChIP-seq pipeline for centromere / repeat enrichment

One pipeline for ChIP-seq experiments whose question lives in repetitive
regions (centromeric alpha-satellite higher-order repeats, human satellites),
where standard MAPQ-filtered workflows delete the signal. Two presets cover
the common designs:

- `--preset pol2`: single-end reads, no spike-in, RPKM tracks.
- `--preset rad21`: paired-end reads, exogenous spike-in chromatin (e.g.
  mouse MEF), spike-scaled tracks.

Three tools, run in this order:

| Tool | Scope | Job |
|---|---|---|
| `bin/chip_repeat_pipeline.sh` | per sample | merge lanes → (downsample) → trim → align → dedup → filter → (bigwig) |
| `bin/spikein_scale_bigwigs.sh` | per experiment | spike factors + spike-scaled bigwigs (spike-in runs only) |
| `bin/hor_enrichment.py` | per experiment | per-chromosome enrichment at HORs/HSat + lollipop plots |

## Requirements

`trimmomatic`, `bowtie2`, `samtools` (≥1.10), `deeptools` (brings `pyBigWig`),
`seqtk` (only for `-d`), `python3` + `matplotlib`, optional `fastqc`.
All on Bioconda (`environment.yml` in this repo includes them). Reference
files: bowtie2 index of the target genome (a T2T-class assembly such as
CHM13v2.0; repeat quantification needs complete centromeres), bowtie2 index
of the spike genome for spike-in runs, `chrom.sizes` of the target genome, and
region BEDs (HOR, optionally HSat2/HSat3, e.g. from the CHM13 censat
annotation) **in target-genome coordinates**.

## The design decisions (why the pipeline is the way it is)

1. **Multi-mappers are kept.** The final filter is `-F 2308`
   (unmapped + secondary + supplementary) with **no MAPQ cut**. HOR alpha
   satellite reads are largely MAPQ0/1; a MAPQ filter deletes the signal you
   are quantifying. Consequence: reads sit at their best (sometimes arbitrary
   within-array) position, so interpret enrichment at the region/chromosome
   level, not base-pair level.
2. **One normalisation, never two.** No spike-in → `RPKM`. Spike-in →
   `--scaleFactor F --normalizeUsing None`. deepTools multiplies
   `--scaleFactor` on top of `--normalizeUsing`, so combining spike factors
   with RPKM double-corrects depth and distorts between-sample comparisons.
   The group script therefore refuses to do it.
3. **Spike factor** = reference / sample spike count, with spike count =
   duplicate-removed primary mapped reads on the spike genome (same trimmed
   reads, same aligner settings as the target). Reference is the lowest spike
   count (Active Motif convention), `1e6` with `-M` (ChIP-Rx), or pinned with
   `-r` (the choice only multiplies all factors by one constant, so
   between-sample comparisons are unchanged).
4. **Enrichment is base-weighted.** Region signal = Σ(signal × bases) read
   directly from the bigwigs (pyBigWig), not `bedmap --sum` over bedgraph
   runs: bigWigToBedGraph run-length-encodes equal-value stretches, and
   summing per-line values under-counts flat tracks (inputs) relative to spiky
   ones (IPs), biasing the ratio.
5. **Alignment**: `bowtie2 --very-sensitive -N 0` (override with `-B`).
6. **Enrichment pairs are yours to choose.** IP vs matched input, IP vs IgG,
   or IP vs IP (e.g. treated vs control). Always ratio two bigwigs with the
   *same* normalisation.

## Worked example 1: single-end, no spike-in (pol2 preset)

RNAP2 ChIP-seq in quiescent cells, TOP2 inhibitor (ICRF-187) vs DMSO, with an
IgG control; four lanes per library:

```bash
A=/path/all_adaptors.fa
IDX=/path/bowtie2_chm13/chm13v2.0

# optional -d: downsample every depth-matched sample to the same read count
bin/chip_repeat_pipeline.sh -P pol2 -n POL2_DMSO -a $A -g $IDX -d 6498399 \
  -1 DMSO_L001_R1.fastq.gz,DMSO_L002_R1.fastq.gz,DMSO_L003_R1.fastq.gz,DMSO_L004_R1.fastq.gz \
  -o results -p 10
bin/chip_repeat_pipeline.sh -P pol2 -n POL2_ICRF -a $A -g $IDX -d 6498399 \
  -1 ICRF_L001_R1.fastq.gz,ICRF_L002_R1.fastq.gz,ICRF_L003_R1.fastq.gz,ICRF_L004_R1.fastq.gz \
  -o results -p 10
bin/chip_repeat_pipeline.sh -P pol2 -n UNT_IgG -a $A -g $IDX \
  -1 IgG_L001_R1.fastq.gz,IgG_L002_R1.fastq.gz,IgG_L003_R1.fastq.gz,IgG_L004_R1.fastq.gz \
  -o results -p 10
```

Lanes of one library are merged before trimming. Then the enrichment:

```bash
cat > pairs_pol2.tsv <<'EOF'
ICRF187_vs_DMSO	results/POL2_ICRF/tracks/POL2_ICRF.rpkm.bw	results/POL2_DMSO/tracks/POL2_DMSO.rpkm.bw
ICRF187_vs_IgG	results/POL2_ICRF/tracks/POL2_ICRF.rpkm.bw	results/UNT_IgG/tracks/UNT_IgG.rpkm.bw
DMSO_vs_IgG	results/POL2_DMSO/tracks/POL2_DMSO.rpkm.bw	results/UNT_IgG/tracks/UNT_IgG.rpkm.bw
EOF

bin/hor_enrichment.py --pairs pairs_pol2.tsv \
  --regions HOR=hor.sorted.bed --regions HSat2=hsat2.bed --regions HSat3=hsat3.bed \
  --chrom-sizes chm13v2.0.chrom.sizes --out results/_enrichment_pol2
```

This yields, per pair: `HOR`, `HSat2`, `HSat3` and `outside_HOR` (chromosome
arm) lollipop panels plus `enrichment.tsv` with the underlying numbers
(including an `ALL` genome-wide row per set).

## Worked example 2: paired-end with spike-in (rad21 preset)

Cohesin (RAD21) ChIP-seq across four conditions, each sample carrying mouse
spike-in chromatin, with matched inputs:

```bash
A=/path/all_adaptors.fa
IDX=/path/target_genome/bowtie2/genome   # in our RAD21 analysis: a CHM13-derived haploid RPE1 assembly
MIDX=/path/mouse/GRCm39                  # spike genome index

for s in IP_cond1 IP_cond2 IP_cond3 IP_cond4 IN_cond1 IN_cond2 IN_cond3 IN_cond4; do
  bin/chip_repeat_pipeline.sh -P rad21 -n $s -a $A -g $IDX -s $MIDX \
    -1 ${s}_R1.fastq.gz -2 ${s}_R2.fastq.gz \
    -o results -p 16
done

# group step: factors + spike-scaled bigwigs (all samples at once)
bin/spikein_scale_bigwigs.sh -o results -p 16
```

`results/_spike/scale_factors.tsv` lists the factors; the bigwigs are built
with `--scaleFactor F --normalizeUsing None`.

```bash
cat > pairs_rad21.tsv <<'EOF'
cond1_vs_input	results/IP_cond1/tracks/IP_cond1.spike.bw	results/IN_cond1/tracks/IN_cond1.spike.bw
cond2_vs_input	results/IP_cond2/tracks/IP_cond2.spike.bw	results/IN_cond2/tracks/IN_cond2.spike.bw
EOF

bin/hor_enrichment.py --pairs pairs_rad21.tsv \
  --regions HOR=hor.sorted.bed \
  --chrom-sizes genome.chrom.sizes --out results/_enrichment_rad21
```

> The region BEDs must match the assembly the reads were aligned to. If you
> align to a custom or patched assembly, confirm its coordinates agree with
> the annotation's source assembly (or lift the annotation over first). For
> example, the CHM13 HOR bed stays valid on a CHM13-derived RPE1 assembly
> only if the assembly preserves CHM13 coordinates (SNV substitutions are
> fine, indels are not).

## Outputs

```
results/<sample>/
  trim/    merged + trimmed FASTQs (deleted intermediates unless -k)
  align/   <sample>.target.final.bam (+.bai)   # dedup, -F 2308, multi-mappers kept
  spike/   spike_count.txt                      # dedup-mapped spike reads
  tracks/  <sample>.rpkm.bw | <sample>.spike.bw
  qc/      fastqc (if installed)
  logs/    trimmomatic, bowtie2 (target+spike), markdup, bamCoverage
  stats.tsv
results/_spike/scale_factors.tsv
results/_enrichment*/enrichment.tsv + <pair>.<set>.png
```

`stats.tsv` fields (labels match the commands):

| field | command |
|---|---|
| `*_total_primary` | `samtools view -c -F 2304` on the raw BAM (all primary records) |
| `*_primary_mapped` | `samtools view -c -F 2308` on the raw BAM (mapped, duplicates still in) |
| `*_dedup_mapped` | `samtools view -c` on the final BAM (mapped, dedup-removed, multi-mappers in) |

## Enrichment details

- Per chromosome and per region set: score = IP signal / control signal, both
  as base-weighted sums over the same intervals, so region length cancels in
  the ratio and scores are comparable across chromosomes with different-sized
  arrays.
- `outside_<set>` is the per-chromosome complement of the first region set
  (chromosome arms), computed from `chrom.sizes`; `chrM` is excluded by
  default (`--exclude-chroms`).
- The per-chromosome enrichment-score and lollipop-plot approach follows
  Saayman et al. 2023 (Mol Cell 83:523–538). `--compiled-pairs` accepts that
  workflow's compiled bedmap files (`chrom start end score`) for
  cross-checking old results; unlike the original script, rows are joined by
  chromosome name rather than row order.

## QC

Group QC (fingerprints, correlation, PCA, MultiQC): the `qc.sh` script from
our companion repo works on the same results layout
(github.com/santoshbiowarrior333/chromatin-pipelines). Note its fingerprint
uses `--minMappingQuality 30`, which is fine as generic QC but by construction
blind to the MAPQ0 repeat signal this pipeline keeps.
