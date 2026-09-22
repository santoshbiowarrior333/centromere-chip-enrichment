# How to run the chip_repeat pipeline, step by step

This is the long version. It walks through every step and explains what each
option does and when you would change it. The short version with worked
examples is in `docs/chip_repeat.md`.

## Before you start

You need five things ready. Check them off first, it saves debugging later.

1. The conda environment:

   ```bash
   conda env create -f environment.yml
   conda activate chromatin
   ```

   On a cluster with modules instead, load: trimmomatic, bowtie2, samtools,
   deeptools, seqtk. Python needs matplotlib and pyBigWig (deeptools brings
   pyBigWig with it).

2. A bowtie2 index of your target genome. Build once if you do not have one:

   ```bash
   bowtie2-build genome.fasta /path/index_prefix
   ```

   Use a T2T style assembly (for human, CHM13v2.0). Older assemblies have
   gaps at centromeres, so there is nothing to map to. A T2T RPE1 genome
   also works, and in general any bowtie2-indexed genome does; keep your
   region BEDs in that assembly's coordinates.

3. If your samples carry spike-in chromatin: a bowtie2 index of the spike
   genome (for mouse MEF spike-in, the mouse genome).

4. A chrom.sizes file for the target genome. Two columns, chromosome name and
   length. `samtools faidx genome.fasta` then `cut -f1,2 genome.fasta.fai`
   gives you one.

5. Your region BED files (HOR at minimum, HSat2 and HSat3 if you want those
   panels) and the adaptor fasta for Trimmomatic. The BED coordinates must
   belong to the same assembly you align to.

## Step 1. Run the per-sample pipeline

One command per sample. IPs, inputs and IgG all go through the same command,
only the FASTQs and the name change.

Single-end sample, no spike-in:

```bash
bin/chip_repeat_pipeline.sh -P pol2 -n POL2_DMSO -a all_adaptors.fa \
  -g /path/chm13_index \
  -1 L001_R1.fastq.gz,L002_R1.fastq.gz,L003_R1.fastq.gz,L004_R1.fastq.gz \
  -o results -p 10
```

Paired-end sample with spike-in:

```bash
bin/chip_repeat_pipeline.sh -P rad21 -n IP_rep1 -a all_adaptors.fa \
  -g /path/target_index -s /path/mouse_index \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -o results -p 16
```

What each option gives you:

| Option | What it does | When to change it |
|---|---|---|
| `-n` | Sample name. Everything for this sample lands in `results/<name>/` and files are named after it. | Always set it. Keep names shell-safe, no spaces. |
| `-1` | R1 FASTQ files, comma separated. If you give several, they are treated as lanes of one library and merged before trimming. | Always set it. |
| `-2` | R2 FASTQ files. Giving this switches the whole run to paired-end (trimming, alignment, fixmate). Leaving it out means single-end. | Set for PE data, omit for SE. |
| `-g` | Bowtie2 index prefix of the target genome. | Always set it. |
| `-s` | Bowtie2 index prefix of the spike genome. Turns on spike-in mode: the same trimmed reads are also aligned to this genome, dedup-removed mapped reads are counted into `spike/spike_count.txt`, and the bigwig step is skipped until Step 2. | Set only for spike-in samples. Use the same `-s` for every sample of the experiment. |
| `-a` | Adaptor fasta for Trimmomatic ILLUMINACLIP. | Always set it. Same file for all samples. |
| `-P` | Preset. `pol2` sets CROP:70 MINLEN:40 (short SE reads). `rad21` sets CROP:72 MINLEN:50 (PE reads). Both keep HEADCROP:15. | Pick the one matching your read layout. Your own `-H`, `-c`, `-m` values always win over the preset. |
| `-H` | Trimmomatic HEADCROP, bases cut from the 5 prime end. Default 15. | Lower it if your library prep does not need the hard 5 prime clip. Check the fastqc per-base content plot. |
| `-c` | Trimmomatic CROP, read length kept after HEADCROP. Default 72. | Match your read length. Reads shorter than CROP pass through unchanged. |
| `-m` | Trimmomatic MINLEN, reads shorter than this after trimming are dropped. Default 50. | Lower for very short reads, raise for long ones. |
| `-d` | Downsample every FASTQ to this many reads with seqtk before trimming, seed 100 so mates stay paired and reruns give the same reads. 0 means off. | Use only when you want depth-matched samples without spike-in, and give the same number to every sample you compare. With spike-in leave it off, the scale factors handle depth. |
| `-B` | Bowtie2 preset. Default `very-sensitive`. | Only change to reproduce an older mapping. Keep it identical across samples you compare. |
| `-z` | Bigwig bin size in bp. Default 50. | Bigger bins give smaller, smoother files. 50 is fine for HOR work. |
| `-o` | Output root. Default `results`. | One root per experiment keeps the group steps simple. |
| `-p` | Threads. Default 8. | Match your allocation. |
| `-k` | Keep intermediate BAMs and merged FASTQs. Off by default to save disk. | Turn on only when debugging a step. |

What the run does, in order: merge lanes, optional downsample, Trimmomatic,
fastqc if installed, bowtie2 to the target genome, fixmate (PE only), sort,
duplicate removal with markdup, filter with `-F 2308`. That filter drops
unmapped, secondary and supplementary reads and nothing else. There is no
MAPQ filter on purpose: centromeric reads are multi-mappers and a MAPQ cut
deletes exactly the signal this pipeline is for.

Then, if there is no `-s`, it writes an RPKM bigwig and you are done with this
sample. If there is a `-s`, it aligns to the spike genome, writes the spike
count and stops. The bigwig for spike-in samples comes from Step 2, because
the scale factor depends on all samples of the experiment.

Check after each run: open `results/<name>/stats.tsv`. It has the read counts
with honest names (total primary reads, mapped with duplicates, mapped after
dedup). The bowtie2 alignment rate is in `logs/bowtie2.target.log`. For a
human ChIP expect an overall rate in the high 90s and a large multi-mapper
fraction, that is normal here.

## Step 2. Spike-in samples only: scale factors and scaled bigwigs

Run once, after every sample of the experiment has finished Step 1:

```bash
bin/spikein_scale_bigwigs.sh -o results -p 16
```

| Option | What it does | When to change it |
|---|---|---|
| `-o` | The same results root as Step 1. The script finds every sample that has a `spike/spike_count.txt`. | Always set it. |
| `-M` | Use 1e6 as the reference (ChIP-Rx style, reads per million spike reads) instead of the lowest sample. | Matter of convention. Between-sample comparisons are identical either way. |
| `-r` | Pin the reference spike count to a number you give. | Use to reproduce factors from an earlier analysis. Also just a constant, comparisons do not change. |
| `-z` | Bigwig bin size. Default 50. | Keep equal to Step 1 `-z`. |
| `-p` | Threads. | As available. |
| trailing names | Restrict to listed samples. | Normally leave empty so all samples get factors from the same reference. |

It writes `results/_spike/scale_factors.tsv` (sample, spike count, reference,
factor) and one `tracks/<sample>.spike.bw` per sample, built with
`--scaleFactor F --normalizeUsing None`.

Do not make spike-in bigwigs by hand with RPKM plus a scale factor. deepTools
multiplies the two, which corrects depth twice and breaks the comparison
between samples. That is the main thing this script exists to prevent.

## Step 3. Enrichment scores and lollipop plots

First write a pairs file. This file is not produced by the pipeline, you
make it yourself in a text editor or with cat, one line per comparison,
pointing at the bigwigs from Steps 1 and 2. Three tab separated columns: a
label for the comparison, the IP bigwig, the control bigwig.

```
ICRF187_vs_DMSO	results/POL2_ICRF/tracks/POL2_ICRF.rpkm.bw	results/POL2_DMSO/tracks/POL2_DMSO.rpkm.bw
```

The control can be a matched input, IgG, or another IP such as the DMSO
condition. The one rule: both bigwigs in a pair must carry the same
normalisation, so rpkm with rpkm, spike with spike.

Then:

```bash
bin/hor_enrichment.py --pairs pairs.tsv \
  --regions HOR=hor.sorted.bed --regions HSat2=hsat2.bed \
  --chrom-sizes genome.chrom.sizes --out results/_enrichment
```

| Option | What it does | When to change it |
|---|---|---|
| `--pairs` | The pairs file above. Every pair is scored against every region set. | Always set it (or use `--compiled-pairs` instead, below). |
| `--regions NAME=BED` | A region set. Repeat the flag for more sets, each becomes its own panel and its own rows in the table. | HOR at minimum. Add HSat2, HSat3 or anything else you want scored. |
| `--chrom-sizes` | Chromosome lengths of the target genome. Needed for the outside set. | Always set it. |
| `--out` | Output folder. | One folder per experiment. |
| `--no-outside` | Skip the outside set. | Only if you do not want the chromosome arm control. |
| `--outside-from` | Which region set to complement for the outside panel. Default is the first `--regions`. | Set when HOR is not your first set. |
| `--exclude-chroms` | Comma list dropped everywhere. Default `chrM`. | Add chrY for female lines, or contigs you do not want plotted. |
| `--colors` | Plot colours, cycled per region set. | Cosmetic. |
| `--compiled-pairs` | Score the old bedmap compiled files (chrom start end score) instead of bigwigs. Rows are matched by chromosome name. | Use to cross-check old results against the new base-weighted ones. |

What the score is: for each chromosome and each region set, the base-weighted
signal sum of the IP inside those intervals divided by the same sum for the
control. A score of 1 means no enrichment, the dashed line in the plots. The
outside set is the complement of the HOR intervals on each chromosome, so it
acts as the arm control: real centromeric enrichment shows scores above 1
within HORs and around 1 outside.

You get `enrichment.tsv` with every number (plus an ALL row per set, the
genome-wide score) and one PNG per pair per set, chromosomes sorted by score.

## Step 4. Sanity checks before you trust the result

1. `scale_factors.tsv`: factors should be within roughly 0.2 to 5. A wild
   factor means a spike alignment problem, check `logs/bowtie2.spike.log`.
2. `enrichment.tsv`: control signal column should not be zero or near zero
   for scored chromosomes. NA in the score column means the control had no
   signal there.
3. The outside panel should sit near 1. If arms are strongly enriched too,
   the story is global, not centromeric.
4. IgG or input pairs (for example DMSO vs IgG) tell you how much apparent
   HOR enrichment exists without any treatment. Look at it once per
   experiment.
5. Group QC: the `qc.sh` script from our companion repo
   (github.com/santoshbiowarrior333/chromatin-pipelines) runs fingerprints,
   correlation, PCA and MultiQC on the same results layout.

## Common problems

- "ERROR: trimmomatic not in PATH": activate the environment or load the
  module before running.
- "ERROR: seqtk needed for -d": seqtk is only required when downsampling.
  Install it or drop `-d`.
- bowtie2 "does not exist or is not a Bowtie 2 index": `-g` and `-s` take the
  index prefix, the path up to but without `.1.bt2`.
- Enrichment table is empty or misses chromosomes: chromosome names differ
  between the BED, the chrom.sizes and the bigwig (chr1 vs 1). Make them
  match, the script only scores names present in the bigwig.
- Spike samples have no bigwig after Step 1: expected. Run Step 2.
- Windows editing note: scripts are LF via .gitattributes. If a script dies
  with `$'\r': command not found`, it picked up CRLF somewhere, run
  `dos2unix` on it or re-checkout.
