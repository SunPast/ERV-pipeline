#!/usr/bin/env python3
"""Validate differences across the three annotations and emit one table (written into the report, to avoid arguments over definitions).

Metrics:
  canonical_loci        -- number of fragments in the unified reference (TEtranscripts rmsk)
  squire_loci           -- number of fragments in the SQuIRE BED
  subfamily_overlap_med -- median per-subfamily locus overlap rate
  subfamily_ge95        -- number of subfamilies with overlap rate >= 95%
  provirus_loci         -- number of loci in the provirus annotation (Telescope 14,968)
  provirus_ltr_share    -- fraction of provirus loci containing LTR elements (prototype/oneside)
"""
import collections
import csv
import sys

map_in, squire_bed, provirus_gtf, out = sys.argv[1:5]

canon = collections.defaultdict(set)
for r in csv.DictReader(open(map_in), delimiter="\t"):
    canon[r["subfamily"]].add((r["chr"], r["start"], r["end"]))

sq = collections.defaultdict(set)
for line in open(squire_bed):
    f = line.split("\t")
    rep = f[3].split("|")[3].split(":")[0]
    if rep in canon:
        sq[rep].add((f[0], f[1], f[2]))

fracs = []
for name, cset in canon.items():
    s = sq.get(name, set())
    fracs.append(len(cset & s) / len(cset) if cset else 0)
fracs.sort()
med = fracs[len(fracs) // 2] if fracs else 0
ge95 = sum(1 for x in fracs if x >= 0.95)

prov = 0
prov_with_ltr = 0
prov_models = set()
ltr_names = set()
with open(provirus_gtf) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        prov += 1
        if 'category "prototype"' in f[8] or 'category "oneside"' in f[8]:
            prov_with_ltr += 1
        i = f[8].find('intModel "')
        if i >= 0:
            prov_models.add(f[8][i + 10: f[8].find('"', i + 10)])

rows = [
    ("canonical_loci", sum(len(v) for v in canon.values())),
    ("canonical_subfamilies", len(canon)),
    ("squire_loci", sum(len(v) for v in sq.values())),
    ("subfamily_overlap_med", round(med, 4)),
    ("subfamily_ge95", ge95),
    ("provirus_loci", prov),
    ("provirus_with_ltr_category", prov_with_ltr),
    ("provirus_internal_models", len(prov_models)),
]
with open(out, "w") as fo:
    fo.write("metric\tvalue\n")
    for k, v in rows:
        fo.write(f"{k}\t{v}\n")
print("\n".join(f"{k}\t{v}" for k, v in rows))
