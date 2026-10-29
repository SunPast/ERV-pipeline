#!/usr/bin/env python3
"""Estimate one size factor per sample for the gene expression tables.

Why a separate file: a raw count is only comparable across samples after
dividing by a sample-level factor, and the ERV size factors must NOT be reused
here. The ERV factors are estimated from ERV counts; applying them to gene
counts would import the ERV composition of each library into the gene layer.

Method: DESeq2's median-of-ratios (Love et al. 2014), estimated across all
samples of all cohorts, using genes with a nonzero count in every sample.
If that set is too small (< --min-genes), the script falls back to library-size
normalisation and records `method = library_size` so the two are never confused.

Usage:
  compute_gene_factors.py OUTDIR SAMPLES_TSV

Output:
  OUTDIR/gene_expression/sample_gene_factor.tsv
    cohort, sample, total_gene_count, n_genes_detected, gene_size_factor, method

The value written is the factor a count must be divided by:
  norm_count = count / gene_size_factor
"""
import csv
import gzip
import os
import sys

import numpy as np
import pandas as pd

outdir, samples_tsv = sys.argv[1:3]

samples = list(csv.DictReader(open(samples_tsv), delimiter="\t"))
order = [(s["cohort"], s["sample"]) for s in samples]

mat_dir = os.path.join(outdir, "gene_expression", "matrices")
blocks = []
for cohort in sorted({c for c, _ in order}):
    path = os.path.join(mat_dir, f"{cohort}.counts.tsv.gz")
    if not os.path.exists(path):
        continue
    # dtype must be given per data column: passing a scalar dtype also applies
    # it to the index, which is a gene id string and would fail to parse
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    df = pd.read_csv(path, sep="\t", index_col=0,
                     dtype={c: np.float32 for c in header[1:]})
    blocks.append(df)

if not blocks:
    sys.exit("no matrices found")

mat = pd.concat(blocks, axis=1)
cols = [(c, s) for c, s in order if s in mat.columns]
mat = mat[[s for _, s in cols]]

total = mat.sum(axis=0)
detected = (mat > 0).sum(axis=0)

all_pos = mat[(mat > 0).all(axis=1)]
method = "deseq2_median_of_ratios"
if all_pos.shape[0] >= 100:
    log_geo = np.log(all_pos).mean(axis=1)
    ratios = np.log(all_pos).sub(log_geo, axis=0)
    factors = np.exp(ratios.median(axis=0))
else:
    method = "library_size"
    factors = total / total.median()

out_path = os.path.join(outdir, "gene_expression", "sample_gene_factor.tsv")
with open(out_path, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["cohort", "sample", "total_gene_count", "n_genes_detected",
                "gene_size_factor", "method"])
    for (cohort, sample) in cols:
        w.writerow([cohort, sample, int(total[sample]), int(detected[sample]),
                    f"{factors[sample]:.6f}", method])

print(f"[factors] {len(cols)} samples, method={method}, "
      f"genes with count>0 in all samples: {all_pos.shape[0]}, -> {out_path}",
      file=sys.stderr)
