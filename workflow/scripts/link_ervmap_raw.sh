#!/bin/bash
# Symlink existing ERVmap results into <OUTDIR>/ervmap_raw/ so that rule
# ervmap_import can find them. ERVmap is NOT rerun: its BED is the method.
#
# Source layout (historical): <RAW_ROOT>/<cohort>/ERVmap/<sample>ERVresults.txt
# Target layout (pipeline):   <OUTDIR>/ervmap_raw/<cohort>/<sample>ERVresults.txt
#
# Usage: link_ervmap_raw.sh OUTDIR [RAW_ROOT] [SAMPLES_TSV]
set -uo pipefail
OUTDIR=${1:?usage: link_ervmap_raw.sh OUTDIR [RAW_ROOT] [SAMPLES_TSV]}
RAW=${2:-/data/home/dingjia/ERV-data/rerun}
SAMP=${3:-config/samples.tsv}

n=0; miss=0; : > "$OUTDIR/ervmap_raw_missing.txt"
while IFS=$'\t' read -r cohort sample _rest; do
  [ "$cohort" = "cohort" ] && continue
  src="$RAW/$cohort/ERVmap/${sample}ERVresults.txt"
  # cohort label != data directory for some cohorts (e.g. EGAD00001003977 lives
  # under EGA_IMvigor210), so fall back to a search by sample name.
  if [ ! -f "$src" ]; then
    cand=$(find "$RAW" -maxdepth 3 -name "${sample}ERVresults.txt" 2>/dev/null | head -1)
    [ -n "$cand" ] && src="$cand"
  fi
  if [ -f "$src" ]; then
    mkdir -p "$OUTDIR/ervmap_raw/$cohort"
    ln -sfn "$src" "$OUTDIR/ervmap_raw/$cohort/${sample}ERVresults.txt"
    n=$((n+1))
  else
    echo -e "$cohort\t$sample" >> "$OUTDIR/ervmap_raw_missing.txt"
    miss=$((miss+1))
  fi
done < "$SAMP"
echo "linked=$n missing=$miss (see $OUTDIR/ervmap_raw_missing.txt)"
