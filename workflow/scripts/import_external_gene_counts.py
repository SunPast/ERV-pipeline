#!/usr/bin/env python3
"""Bring an externally computed gene-count matrix into the pipeline.

Some cohorts may already have per-cohort count matrices produced outside this
repository. Re-running featureCounts on them is wasted compute, but silently
mixing numbers of unknown provenance into the database is worse. This script
makes the reuse explicit: it converts those matrices into the pipeline's own
long-table format and records where every value came from.

Input
    --dir       directory of wide matrices, one per cohort
    --pattern   filename pattern, must contain '{cohort}'
    --samples   the pipeline sample sheet (cohort, sample)
    --out       output long table (gzipped TSV)

Output columns
    cohort, sample, gene_id, count, source

Only count > 0 rows are written; the database keeps counts sparse, and the
number of all-zero genes is large.

Samples present in the matrices but absent from the sample sheet are ignored
(they belong to the other project), and samples in the sheet without a matrix
column are reported so the caller knows what still has to be recomputed.
"""
import argparse
import csv
import gzip
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--pattern", default="{cohort}_merged_counts.csv")
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--source-name", default="external")
    args = ap.parse_args()

    with open(args.samples) as fh:
        sheet = {r["sample"]: r["cohort"] for r in csv.DictReader(fh, delimiter="\t")}
    wanted = set(sheet)

    written = 0
    seen = set()
    for fname in sorted(os.listdir(args.dir)):
        if not fname.endswith(".csv"):
            continue
        path = os.path.join(args.dir, fname)
        with open(path) as fh:
            header = next(csv.reader(fh))
            cols = header[1:]
            keep = [(i, c) for i, c in enumerate(cols, start=1) if c in wanted]
            if not keep:
                continue
            n_before = written
            with gzip.open(args.out, "at") as out:
                if out.tell() == 0:
                    out.write("cohort\tsample\tgene_id\tcount\tsource\n")
                for row in csv.reader(fh):
                    gene = row[0]
                    for i, sample in keep:
                        try:
                            v = int(float(row[i] or 0))
                        except (ValueError, IndexError):
                            continue
                        if v > 0:
                            out.write(f"{sheet[sample]}\t{sample}\t{gene}\t{v}\t{args.source_name}\n")
                            written += 1
            seen.update(c for _, c in keep)
            print(f"{fname}: {len(keep)} samples, {written - n_before} rows", file=sys.stderr)

    missing = sorted(wanted - seen)
    print(f"total rows={written} samples_from_external={len(seen)} "
          f"samples_not_covered={len(missing)}", file=sys.stderr)
    if missing:
        with open(args.out + ".missing.tsv", "w") as fh:
            for s in missing:
                fh.write(f"{sheet[s]}\t{s}\n")


if __name__ == "__main__":
    main()
