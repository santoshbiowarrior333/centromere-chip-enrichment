#!/usr/bin/env bash
#
# spikein_scale_bigwigs.sh - group step: spike-in scale factors + scaled bigwigs
#
# Run AFTER chip_repeat_pipeline.sh has finished for ALL samples of an
# experiment (IPs and inputs). It:
#   1. collects each sample's spike-in dedup-mapped read count
#      (results/<sample>/spike/spike_count.txt)
#   2. computes one scale factor per sample:
#        default : factor = min(spike counts) / spike_count   (Active Motif /
#                  Orlando; the reference sample gets factor 1, all others <1)
#        -M      : factor = 1e6 / spike_count                 (ChIP-Rx, RRPM)
#   3. writes results/_spike/scale_factors.tsv
#   4. builds tracks/<sample>.spike.bw with
#        bamCoverage --scaleFactor F --normalizeUsing None
#      (never RPKM here: the spike factor already carries the depth correction;
#       combining the two double-corrects depth and distorts comparisons)
#
# Usage:
#   spikein_scale_bigwigs.sh -o results [-M | -r REF_COUNT] [-z 50] [-p 8] [sample ...]
#
# Options:
#   -o  results root written by chip_repeat_pipeline.sh (default results)
#   -M  per-million reference (ChIP-Rx) instead of min-sample reference
#   -r  explicit reference spike count, e.g. to reproduce the factors of a
#       previous analysis. Changing the reference multiplies every factor by
#       the same constant, so between-sample comparisons are identical.
#   -z  bigwig bin size (default 50)
#   -p  threads (default 8)
#   positional: restrict to these sample names (default: every sample with a
#   spike count under -o)
#
set -euo pipefail

OUTROOT="results" ; PERMILLION=0 ; REF_OVERRIDE="" ; BINSIZE=50 ; THREADS=8

usage() { grep '^#' "$0" | sed 's/^# \{0,1\}//' | sed -n '2,36p'; exit 1; }

while getopts "o:Mr:z:p:h" opt; do
  case $opt in
    o) OUTROOT=$OPTARG ;;
    M) PERMILLION=1 ;;
    r) REF_OVERRIDE=$OPTARG ;;
    z) BINSIZE=$OPTARG ;;
    p) THREADS=$OPTARG ;;
    h|*) usage ;;
  esac
done
shift $((OPTIND-1))

command -v bamCoverage >/dev/null 2>&1 || { echo "ERROR: bamCoverage not in PATH"; exit 1; }
[[ -d "$OUTROOT" ]] || { echo "ERROR: $OUTROOT not found"; exit 1; }

declare -a SAMPLES COUNTS BAMS
if [[ $# -gt 0 ]]; then
  want=("$@")
else
  want=()
  for d in "$OUTROOT"/*/; do
    [[ -f "$d/spike/spike_count.txt" ]] && want+=("$(basename "$d")")
  done
fi
[[ ${#want[@]} -ge 1 ]] || { echo "ERROR: no samples with spike/spike_count.txt under $OUTROOT"; exit 1; }

for s in "${want[@]}"; do
  cfile="$OUTROOT/$s/spike/spike_count.txt"
  bam="$OUTROOT/$s/align/$s.target.final.bam"
  [[ -f "$cfile" ]] || { echo "ERROR: $cfile missing"; exit 1; }
  [[ -f "$bam"   ]] || { echo "ERROR: $bam missing"; exit 1; }
  c=$(<"$cfile")
  [[ "$c" =~ ^[0-9]+$ && "$c" -gt 0 ]] || { echo "ERROR: bad spike count for $s: '$c'"; exit 1; }
  SAMPLES+=("$s"); COUNTS+=("$c"); BAMS+=("$bam")
done

# reference
REF=""
if [[ -n "$REF_OVERRIDE" ]]; then
  [[ "$REF_OVERRIDE" =~ ^[0-9]+$ ]] || { echo "ERROR: -r must be an integer"; exit 1; }
  REF=$REF_OVERRIDE
  echo "reference: user-supplied = $REF"
elif [[ $PERMILLION -eq 1 ]]; then
  REF=1000000
  echo "reference: 1e6 spike reads (ChIP-Rx / RRPM)"
else
  for c in "${COUNTS[@]}"; do
    [[ -z "$REF" || "$c" -lt "$REF" ]] && REF=$c
  done
  echo "reference: lowest spike count = $REF"
fi

mkdir -p "$OUTROOT/_spike"
FACTORS="$OUTROOT/_spike/scale_factors.tsv"
{
  echo -e "sample\tspike_dedup_mapped\treference\tscale_factor"
  for i in "${!SAMPLES[@]}"; do
    f=$(awk -v r="$REF" -v c="${COUNTS[$i]}" 'BEGIN{printf "%.10f", r/c}')
    echo -e "${SAMPLES[$i]}\t${COUNTS[$i]}\t${REF}\t${f}"
  done
} > "$FACTORS"
echo "wrote $FACTORS"
column -t "$FACTORS" 2>/dev/null || cat "$FACTORS"

for i in "${!SAMPLES[@]}"; do
  s=${SAMPLES[$i]}
  f=$(awk -F'\t' -v s="$s" '$1==s{print $4}' "$FACTORS")
  out="$OUTROOT/$s/tracks/$s.spike.bw"
  echo "[$(date +%T)] bamCoverage $s (scaleFactor $f, normalizeUsing None)"
  bamCoverage -b "${BAMS[$i]}" -o "$out" \
    --scaleFactor "$f" --normalizeUsing None \
    --binSize "$BINSIZE" -p "$THREADS" \
    2> "$OUTROOT/$s/logs/bamCoverage.spike.log"
done

echo "done. Spike-scaled bigwigs: $OUTROOT/<sample>/tracks/<sample>.spike.bw"
