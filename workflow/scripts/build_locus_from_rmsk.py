#!/usr/bin/env python3
"""Build locus-level GTF / map tables directly from the UCSC RepeatMasker table.

This replaces the TEtranscripts-GTF-based build: coordinates, subfamily,
family, class, strand AND milliDiv all come from ONE file, so there is no
cross-version mixing. The ERV subfamily set is still defined by the dictionary.

Usage:
  build_locus_from_rmsk.py RMSK_TABLE DICT GTF_OUT MAP_OUT EXON_OUT [CHROMS]

  RMSK_TABLE : UCSC hg38_rmsk.txt (17 columns, 1-based file, 0-based coords)
  DICT       : subfamily dictionary CSV; "-" disables the subfamily filter
  GTF_OUT / EXON_OUT : "-" to skip writing that file
  CHROMS     : regex of allowed chromosomes, default "^chr([0-9]+|X|Y)$"

RMSK columns (0-based): 0 bin, 1 swScore, 2 milliDiv, 3 milliDel, 4 milliIns,
5 genoName, 6 genoStart, 7 genoEnd, 8 genoLeft, 9 strand, 10 repName,
11 repClass, 12 repFamily.

Output map columns:
  locus_id, subfamily, family, rep_class, region_type, chr, start, end,
  strand, milli_div
"""
import csv
import re
import sys

rmsk_in, dict_in, gtf_out, map_out, exon_out = sys.argv[1:6]
chrom_re = re.compile(sys.argv[6] if len(sys.argv) > 6 else r"^chr([0-9]+|X|Y)$")

if dict_in == "-":
    names = None
else:
    names = {r["subfamily"] for r in csv.DictReader(open(dict_in))}

recs = {}
with open(rmsk_in) as fh:
    for line in fh:
        p = line.split()
        if len(p) < 13:
            continue
        chrom, start, end, strand = p[5], int(p[6]), int(p[7]), p[9]
        rep, cls, fam, md = p[10], p[11], p[12], p[2]
        if not chrom_re.match(chrom):
            continue
        if names is not None and rep not in names:
            continue
        key = (rep, chrom, start, end, strand)
        if key in recs:
            continue
        recs[key] = (fam, cls, md)

fo = open(gtf_out, "w") if gtf_out != "-" else None
fe = open(exon_out, "w") if exon_out != "-" else None
with open(map_out, "w") as fm:
    fm.write("locus_id\tsubfamily\tfamily\trep_class\tregion_type\tchr\tstart\tend\tstrand\tmilli_div\n")
    for (rep, chrom, start, end, strand), (fam, cls, md) in sorted(
        recs.items(), key=lambda kv: (kv[0][1], kv[0][2], kv[0][3])
    ):
        locus = f"{rep}__{chrom}_{start}_{end}_{strand}"
        reg = "internal" if rep.endswith("-int") else "ltr"
        attr = (
            f'gene_id "{locus}"; transcript_id "{locus}"; locus "{locus}"; '
            f'repName "{rep}"; repFamily "{fam}"; repClass "{cls}"; '
            f'family_id "{fam}"; class_id "{cls}"; '
            f'subfamily "{rep}"; region_type "{reg}"; milliDiv "{md}";'
        )
        line = "\t".join([chrom, "rmsk", "gene", str(start + 1), str(end),
                          ".", strand, ".", attr]) + "\n"
        if fo:
            fo.write(line)
        if fe:
            fe.write(line.replace("\trmsk\tgene\t", "\trmsk\texon\t"))
        fm.write(f"{locus}\t{rep}\t{fam}\t{cls}\t{reg}\t{chrom}\t{start}\t{end}\t{strand}\t{md}\n")

if fo:
    fo.close()
if fe:
    fe.close()

print(f"loci={len(recs)}", file=sys.stderr)
