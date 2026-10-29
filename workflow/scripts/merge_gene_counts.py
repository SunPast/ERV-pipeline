#!/usr/bin/env python3
"""Merge per-sample featureCounts outputs into the gene expression tables.

Usage:
  merge_gene_counts.py OUTDIR SAMPLES_TSV GTF

Outputs, both derived from the per-sample `_featureCounts.txt` files:

  OUTDIR/gene_expression/matrices/<cohort>.counts.tsv.gz
      wide matrix, one row per gene (Geneid with version), one column per
      sample of that cohort. Kept unfiltered so downstream analyses can apply
      their own rule.

  OUTDIR/gene_expression/gene_expression_standard.csv.gz
      long table for the database loader:
      cohort,sample,gene_id,gene_symbol,length,count
      only rows with count > 0, matching the convention of the ERV standard
      tables (zero-count features are not shipped).

featureCounts writes one header line starting with '#' plus a column header;
the last column is this sample's count because each run has exactly one BAM.
"""
import csv
import gzip
import os
import re
import sys

outdir, samples_tsv, gtf_path = sys.argv[1:4]

# ---- gene metadata from the annotation ----
gene_symbol = {}
gene_len = {}
pat = re.compile(r'gene_id "([^"]+)".*gene_name "([^"]+)"')
for line in open(gtf_path):
    if "\tgene\t" not in line:
        continue
    m = pat.search(line)
    if not m:
        continue
    f = line.rstrip("\n").split("\t")
    gene_symbol[m.group(1)] = m.group(2)
    gene_len[m.group(1)] = int(f[4]) - int(f[3]) + 1

samples = list(csv.DictReader(open(samples_tsv), delimiter="\t"))
by_cohort = {}
for s in samples:
    by_cohort.setdefault(s["cohort"], []).append(s["sample"])

mat_dir = os.path.join(outdir, "gene_expression", "matrices")
os.makedirs(mat_dir, exist_ok=True)
os.makedirs(os.path.join(outdir, "gene_expression"), exist_ok=True)

std_path = os.path.join(outdir, "gene_expression", "gene_expression_standard.csv.gz")
std = gzip.open(std_path, "wt", newline="")
w = csv.writer(std)
w.writerow(["cohort", "sample", "gene_id", "gene_symbol", "length", "count"])

for cohort, smps in sorted(by_cohort.items()):
    cols = {}
    for smp in smps:
        path = os.path.join(outdir, "gene_counts", cohort, smp,
                            f"{smp}_featureCounts.txt")
        if not os.path.exists(path):
            continue
        counts = {}
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if p[0] == "Geneid":
                    continue
                counts[p[0]] = int(float(p[-1]))
        cols[smp] = counts
        for gid, c in counts.items():
            if c > 0:
                w.writerow([cohort, smp, gid, gene_symbol.get(gid, ""),
                            gene_len.get(gid, ""), c])

    if not cols:
        continue
    genes = sorted(set().union(*cols.values()))
    with gzip.open(os.path.join(mat_dir, f"{cohort}.counts.tsv.gz"), "wt") as mh:
        mh.write("\t".join(["Geneid"] + list(cols.keys())) + "\n")
        for gid in genes:
            mh.write("\t".join([gid] + [str(c.get(gid, 0)) for c in cols.values()]) + "\n")

std.close()
print(f"gene expression written to {outdir}/gene_expression", file=sys.stderr)
