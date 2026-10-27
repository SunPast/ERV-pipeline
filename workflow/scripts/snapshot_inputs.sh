#!/bin/bash
# Snapshot the inputs that the SQL side uses to compute norm_count
# (sample -> pseudonym map, and size factors) together with the run products.
# Without this, a regenerated size-factor file silently changes every
# norm_count in the database with no trace of which version was used.
#
# Usage: snapshot_inputs.sh OUTDIR FILE [FILE ...]
set -uo pipefail
OUT=${1:?usage: snapshot_inputs.sh OUTDIR FILE [FILE ...]}
shift
mkdir -p "$OUT"
for f in "$@"; do
  [ -f "$f" ] || { echo "missing: $f" >&2; continue; }
  cp -f "$f" "$OUT/"
done
( cd "$OUT" && md5sum ./* > MD5SUMS )
echo "snapshot written: $OUT"
