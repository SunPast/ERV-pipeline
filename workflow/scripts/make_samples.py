#!/usr/bin/env python3
"""Build config/samples.tsv and record how every sample was resolved.

Sample -> BAM has two sources:
  1. a historical manifest (sample<TAB>bam), produced by the old Telescope run;
  2. a filesystem search by sample name, for samples missing from that manifest.

The second source is not deterministic enough to leave unrecorded, so a
resolution log is always written:
  cohort, sample, method(manifest|glob|glob_multi|missing), n_candidates, bam

Usage:
  make_samples.py OUT_TSV OUT_LOG COHORT_MAP BAM_MANIFEST BAM_ROOT
"""
import csv
import glob
import os
import sys

out_tsv, out_log, cohort_map, bam_manifest, bam_root = sys.argv[1:6]

bams = {}
with open(bam_manifest) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            bams[p[0]] = p[1]

rows = list(csv.DictReader(open(cohort_map)))

seen = set()
n_manifest = n_glob = n_multi = n_missing = 0
with open(out_tsv, "w") as ft, open(out_log, "w") as fl:
    ft.write("cohort\tsample\trun\tbam\tfq1\tfq2\n")
    fl.write("cohort\tsample\tmethod\tn_candidates\tbam\n")
    for r in rows:
        cohort, sample = r["cohort"], r["sample"]
        if sample in seen:
            continue
        seen.add(sample)

        method, cand = "manifest", []
        bam = bams.get(sample, "")
        # A manifest entry that exists but is 0 bytes is a failed STAR run, not a
        # usable BAM (PHS002176 has 4 of them under ERV.map/results). Fall back
        # to the filesystem search so the real alignment is picked up.
        if bam and os.path.exists(bam) and os.path.getsize(bam) > 0:
            n_manifest += 1
        else:
            cand = sorted(glob.glob(
                f"{bam_root}/*/aligned/{sample}/*_Aligned.sortedByCoord.out.bam"))
            if len(cand) == 1:
                bam, method = cand[0], "glob"
                n_glob += 1
            elif len(cand) > 1:
                bam, method = cand[0], "glob_multi"
                n_multi += 1
            else:
                bam, method = "", "missing"
                n_missing += 1

        fl.write(f"{cohort}\t{sample}\t{method}\t{len(cand) if cand else (1 if method == 'manifest' else 0)}\t{bam}\n")
        if bam:
            ft.write(f"{cohort}\t{sample}\t{sample}\t{bam}\t\t\n")

print(f"samples={len(seen)} manifest={n_manifest} glob={n_glob} "
      f"glob_multi={n_multi} missing={n_missing}", file=sys.stderr)
