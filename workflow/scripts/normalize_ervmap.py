#!/usr/bin/env python3
"""Normalize raw ERVmap output into the unified long table (locus level).

Input: <sample>ERVresults.txt -- chr, start, end, ervmap_id, score, strand, count
      locus_catalog.csv      -- platform final locus table (the ERV layer carries the real strand and subfamily)
Output: <sample>.locus.tsv      -- locus_id, chr, start, end, strand, subfamily, family, count
"""
import csv
import sys

raw_in, catalog_in, out = sys.argv[1:4]

cat = {}
with open(catalog_in) as fh:
    for r in csv.DictReader(fh):
        if r["layer"] != "ERV":
            continue
        cat[r["ervmap_id"]] = r

with open(raw_in) as fh, open(out, "w") as fo:
    fo.write("locus_id\tchr\tstart\tend\tstrand\tsubfamily\tfamily\tcount\n")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 7:
            continue
        chrom, start, end, eid, _score, _strand, count = f[:7]
        c = cat.get(eid, {})
        fo.write("\t".join([
            eid, chrom, start, end,
            c.get("strand", "."), c.get("subfamily", ""), c.get("family", ""),
            count,
        ]) + "\n")
