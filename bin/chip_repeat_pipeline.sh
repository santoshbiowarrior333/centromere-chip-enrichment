#!/usr/bin/env bash
#
# chip_repeat_pipeline.sh - unified ChIP-seq pipeline for centromere/repeat analysis
#
# One script for both protocols:
#   POL2-style : single-end, no spike-in, RPKM tracks          (--preset pol2)
#   RAD21-style: paired-end, mouse spike-in, spike-scaled tracks (--preset rad21)
#
# Backbone (identical for both): lane merge -> [seqtk downsample] -> Trimmomatic
# -> bowtie2 --very-sensitive -N 0 -> [fixmate (PE only)] -> sort -> markdup -r
# -> -F 2308 filter (multi-mappers KEPT: no MAPQ filter, required for HOR signal)
# -> bigwig.
#
# Normalisation rules (this is deliberate, do not mix them):
#   no spike-in  -> bamCoverage --normalizeUsing RPKM
#   spike-in     -> bigwig deferred; run spikein_scale_bigwigs.sh after ALL
#                   samples, which uses --scaleFactor F --normalizeUsing None.
#                   NEVER combine --scaleFactor with RPKM (double depth-correction).
#
# Usage:
#   chip_repeat_pipeline.sh -n NAME -g BT2_INDEX -a adaptors.fa \
#       -1 R1_L001.fq.gz,R1_L002.fq.gz [-2 R2_L001.fq.gz,R2_L002.fq.gz] \
#       [-s SPIKE_BT2_INDEX] [-P pol2|rad21] [-d N] [-o results] [-p 8]
#
# Options:
#   -n  sample name (required)
#   -1  comma-separated R1 FASTQ(.gz) list; lanes are merged (required)
#   -2  comma-separated R2 FASTQ(.gz) list -> paired-end mode
#   -g  bowtie2 index prefix, target genome (required)
#   -s  bowtie2 index prefix, spike-in genome (enables spike counting)
#   -a  adaptor fasta for Trimmomatic ILLUMINACLIP (required)
#   -P  preset: pol2 (SE, CROP 70, MINLEN 40) | rad21 (PE, CROP 72, MINLEN 50)
#   -H  Trimmomatic HEADCROP   (default 15)
#   -c  Trimmomatic CROP       (default 72)
#   -m  Trimmomatic MINLEN     (default 50)
#   -d  downsample to N reads with seqtk (seed 100) before trimming; 0=off (default 0)
#   -B  bowtie2 preset (default very-sensitive)
#   -z  bigwig bin size (default 50)
#   -o  output root (default results)
#   -p  threads (default 8)
#   -k  keep intermediate BAM/FASTQ files
#   -h  this help
#
# Output: OUTROOT/NAME/{trim,align,spike,tracks,qc,logs} + OUTROOT/NAME/stats.tsv
#
set -euo pipefail

NAME="" ; R1="" ; R2="" ; IDX="" ; SPIKE_IDX="" ; ADAPTORS=""
PRESET="" ; HEADCROP=15 ; CROP=72 ; MINLEN=50 ; DOWNSAMPLE=0
BT2_PRESET="very-sensitive" ; BINSIZE=50 ; OUTROOT="results" ; THREADS=8 ; KEEP=0

usage() { grep '^#' "$0" | sed 's/^# \{0,1\}//' | sed -n '2,50p'; exit 1; }

while getopts "n:1:2:g:s:a:P:H:c:m:d:B:z:o:p:kh" opt; do
  case $opt in
    n) NAME=$OPTARG ;;
    1) R1=$OPTARG ;;
    2) R2=$OPTARG ;;
    g) IDX=$OPTARG ;;
    s) SPIKE_IDX=$OPTARG ;;
    a) ADAPTORS=$OPTARG ;;
    P) PRESET=$OPTARG ;;
    H) HEADCROP=$OPTARG ;;
    c) CROP=$OPTARG ; CROP_SET=1 ;;
    m) MINLEN=$OPTARG ; MINLEN_SET=1 ;;
    d) DOWNSAMPLE=$OPTARG ;;
    B) BT2_PRESET=$OPTARG ;;
    z) BINSIZE=$OPTARG ;;
    o) OUTROOT=$OPTARG ;;
    p) THREADS=$OPTARG ;;
    k) KEEP=1 ;;
    h|*) usage ;;
  esac
done

# presets only fill defaults the user did not override on the command line
case "$PRESET" in
  pol2)  [[ -z "${CROP_SET:-}"   ]] && CROP=70
         [[ -z "${MINLEN_SET:-}" ]] && MINLEN=40 ;;
  rad21) [[ -z "${CROP_SET:-}"   ]] && CROP=72
         [[ -z "${MINLEN_SET:-}" ]] && MINLEN=50 ;;
  "") : ;;
  *) echo "ERROR: unknown preset '$PRESET' (pol2|rad21)"; exit 1 ;;
esac

[[ -n "$NAME" && -n "$R1" && -n "$IDX" && -n "$ADAPTORS" ]] || usage
[[ -f "$ADAPTORS" ]] || { echo "ERROR: adaptor fasta not found: $ADAPTORS"; exit 1; }

for tool in trimmomatic bowtie2 samtools bamCoverage; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not in PATH"; exit 1; }
done
if [[ "$DOWNSAMPLE" -gt 0 ]]; then
  command -v seqtk >/dev/null 2>&1 || { echo "ERROR: seqtk needed for -d"; exit 1; }
fi

PE=0 ; [[ -n "$R2" ]] && PE=1
S="$OUTROOT/$NAME"
mkdir -p "$S"/{trim,align,spike,tracks,qc,logs}
log() { echo "[$(date +%T)] $*"; }

# ------------------------------------------------------------------ lane merge
log "merging lanes"
merge_lanes() {  # $1 comma list -> $2 out file (gz-safe: gz members concatenate)
  local IFS=','; local files=($1); local f
  for f in "${files[@]}"; do [[ -f "$f" ]] || { echo "ERROR: missing $f"; exit 1; }; done
  cat "${files[@]}" > "$2"
}
R1M="$S/trim/${NAME}_R1.merged.fastq.gz"
case "$R1" in *.gz*) : ;; *) R1M="${R1M%.gz}" ;; esac
merge_lanes "$R1" "$R1M"
if [[ $PE -eq 1 ]]; then
  R2M="$S/trim/${NAME}_R2.merged.fastq.gz"
  case "$R2" in *.gz*) : ;; *) R2M="${R2M%.gz}" ;; esac
  merge_lanes "$R2" "$R2M"
fi

# ------------------------------------------------------------- downsample (opt)
if [[ "$DOWNSAMPLE" -gt 0 ]]; then
  log "downsampling to $DOWNSAMPLE reads (seqtk, seed 100 - same seed keeps mates paired)"
  seqtk sample -s100 "$R1M" "$DOWNSAMPLE" | gzip > "$S/trim/${NAME}_R1.ds.fastq.gz"
  R1M="$S/trim/${NAME}_R1.ds.fastq.gz"
  if [[ $PE -eq 1 ]]; then
    seqtk sample -s100 "$R2M" "$DOWNSAMPLE" | gzip > "$S/trim/${NAME}_R2.ds.fastq.gz"
    R2M="$S/trim/${NAME}_R2.ds.fastq.gz"
  fi
fi

# -------------------------------------------------------------------- trimming
TSTEPS="ILLUMINACLIP:${ADAPTORS}:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 HEADCROP:${HEADCROP} CROP:${CROP} MINLEN:${MINLEN}"
T1="$S/trim/${NAME}_R1.trimmed.fastq.gz"
if [[ $PE -eq 1 ]]; then
  T2="$S/trim/${NAME}_R2.trimmed.fastq.gz"
  log "trimmomatic PE ($TSTEPS)"
  trimmomatic PE -threads "$THREADS" -phred33 "$R1M" "$R2M" \
    "$T1" "$S/trim/${NAME}_R1.unpaired.fastq.gz" \
    "$T2" "$S/trim/${NAME}_R2.unpaired.fastq.gz" \
    $TSTEPS 2> "$S/logs/trimmomatic.log"
else
  log "trimmomatic SE ($TSTEPS)"
  trimmomatic SE -threads "$THREADS" -phred33 "$R1M" "$T1" \
    $TSTEPS 2> "$S/logs/trimmomatic.log"
fi
if command -v fastqc >/dev/null 2>&1; then
  log "fastqc (post-trim)"
  if [[ $PE -eq 1 ]]; then fastqc -q -o "$S/qc" "$T1" "$T2" || true
  else fastqc -q -o "$S/qc" "$T1" || true; fi
fi

# ------------------------------------------------------- align / dedup / filter
align_dedup_filter() {  # $1 bt2 index, $2 outdir, $3 tag -> writes $2/$NAME.$3.final.bam
  local idx=$1 dir=$2 tag=$3
  local raw="$dir/$NAME.$tag.raw.bam"
  local fixm="$dir/$NAME.$tag.fixmate.bam"
  local sorted="$dir/$NAME.$tag.sorted.bam"
  local nodup="$dir/$NAME.$tag.nodup.bam"
  local final="$dir/$NAME.$tag.final.bam"

  log "bowtie2 --$BT2_PRESET -N 0 -> $tag genome"
  if [[ $PE -eq 1 ]]; then
    bowtie2 --"$BT2_PRESET" -N 0 -p "$THREADS" -x "$idx" -q \
      -1 "$T1" -2 "$T2" 2> "$S/logs/bowtie2.$tag.log" \
      | samtools view -b -o "$raw" -
  else
    bowtie2 --"$BT2_PRESET" -N 0 -p "$THREADS" -x "$idx" -q \
      -U "$T1" 2> "$S/logs/bowtie2.$tag.log" \
      | samtools view -b -o "$raw" -
  fi

  if [[ $PE -eq 1 ]]; then
    log "fixmate -m ($tag)"                       # bowtie2 output is name-collated
    samtools fixmate -m -@ "$THREADS" "$raw" "$fixm"
  else
    fixm="$raw"                                   # fixmate is a no-op on SE data
  fi

  log "sort + markdup -r ($tag)"
  samtools sort -@ "$THREADS" -o "$sorted" "$fixm"
  samtools markdup -r -@ "$THREADS" -s "$sorted" "$nodup" 2> "$S/logs/markdup.$tag.log"

  log "filter -F 2308 ($tag; multi-mappers kept)"
  samtools view -b -@ "$THREADS" -F 2308 -o "$final" "$nodup"
  samtools index "$final"

  # honest read accounting (see docs: labels match the commands)
  local total primary_mapped dedup_mapped
  total=$(samtools view -c -F 2304 -@ "$THREADS" "$raw")          # primary records
  primary_mapped=$(samtools view -c -F 2308 -@ "$THREADS" "$raw") # mapped primary, dups included
  dedup_mapped=$(samtools view -c -@ "$THREADS" "$final")          # mapped primary, dups removed
  echo -e "${tag}_total_primary\t${total}"        >> "$S/stats.raw"
  echo -e "${tag}_primary_mapped\t${primary_mapped}" >> "$S/stats.raw"
  echo -e "${tag}_dedup_mapped\t${dedup_mapped}"  >> "$S/stats.raw"

  if [[ $KEEP -eq 0 ]]; then
    rm -f "$raw" "$sorted" "$nodup"
    [[ $PE -eq 1 ]] && rm -f "$fixm"
  fi
}

: > "$S/stats.raw"
align_dedup_filter "$IDX" "$S/align" "target"

SPIKE_COUNT=""
if [[ -n "$SPIKE_IDX" ]]; then
  align_dedup_filter "$SPIKE_IDX" "$S/spike" "spike"
  SPIKE_COUNT=$(awk -F'\t' '$1=="spike_dedup_mapped"{print $2}' "$S/stats.raw")
  echo "$SPIKE_COUNT" > "$S/spike/spike_count.txt"
  if [[ $KEEP -eq 0 ]]; then rm -f "$S/spike/$NAME.spike.final.bam" "$S/spike/$NAME.spike.final.bam.bai"; fi
fi

# ---------------------------------------------------------------------- bigwig
FINAL="$S/align/$NAME.target.final.bam"
if [[ -z "$SPIKE_IDX" ]]; then
  log "bamCoverage RPKM (no spike-in)"
  bamCoverage -b "$FINAL" -o "$S/tracks/${NAME}.rpkm.bw" \
    --normalizeUsing RPKM --binSize "$BINSIZE" -p "$THREADS" \
    2> "$S/logs/bamCoverage.log"
else
  log "spike-in sample: bigwig DEFERRED."
  log "  after all samples finish, run: spikein_scale_bigwigs.sh -o $OUTROOT"
fi

# ----------------------------------------------------------------------- stats
{
  echo -e "sample\tlayout\tfield\tvalue"
  while IFS=$'\t' read -r k v; do
    echo -e "$NAME\t$([[ $PE -eq 1 ]] && echo PE || echo SE)\t$k\t$v"
  done < "$S/stats.raw"
} > "$S/stats.tsv"
rm -f "$S/stats.raw"

log "done: $S"
log "  final BAM : $FINAL"
[[ -z "$SPIKE_IDX" ]] && log "  bigwig    : $S/tracks/${NAME}.rpkm.bw"
[[ -n "$SPIKE_IDX" ]] && log "  spike-in dedup-mapped reads: $SPIKE_COUNT"
log "  stats     : $S/stats.tsv"
