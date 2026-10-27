#!/usr/bin/env python3
"""ERVmap provirus locus -> composition mapping (internal region / 5'LTR / 3'LTR)

Principles (see docs/annotation.md):
  * Primary assignment = the provirus locus's intModel (internal region);
  * LTR names are taken directly from the rmsk fragments inside the locus (no guessing GTF prefixes);
  * No "coverage-fraction threshold" is used to decide primary/secondary: internal region takes priority, LTRs are always extra fields;
  * When an ERVmap locus does not match a provirus locus, fall back to the rmsk fragment (assign the one with the largest overlap).

Input:
  locus_catalog.csv        ERV layer, 3,220 rows (with real strand)
  provirus_gtf             provirus annotation (14,968 loci, with intModel/category/span)
  te_locus_map.tsv         unified rmsk fragment table (subfamily / region_type / coordinates)
Output:
  ervmap_provirus_map.tsv
    ervmap_id, chr, start, end, strand,
    provirus_locus, provirus_category, int_model,
    ltr5_model, ltr3_model, ltr_models, n_ltr,
    fallback_subfamily, source
"""
import bisect
import collections
import csv
import re
import sys

catalog_in, provirus_gtf, locus_map, out = sys.argv[1:5]

# ---- provirus loci ----
prov = collections.defaultdict(list)
with open(provirus_gtf) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene" or "gene_id" not in f[8]:
            continue

        def grab(key):
            i = f[8].find(key + ' "')
            return f[8][i + len(key) + 2: f[8].find('"', i + len(key) + 2)] if i >= 0 else ""

        prov[f[0]].append((int(f[3]) - 1, int(f[4]), f[6], grab("gene_id"),
                           grab("intModel"), grab("category")))
for k in prov:
    prov[k].sort()

# ---- rmsk fragments ----
frag = collections.defaultdict(list)
with open(locus_map) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        frag[r["chr"]].append((int(r["start"]), int(r["end"]), r["subfamily"],
                               r["region_type"], r["strand"]))
for k in frag:
    frag[k].sort()


def hits(lst, s, e):
    i = bisect.bisect_left(lst, (s,)) - 1
    while i < len(lst):
        a, b = lst[i][0], lst[i][1]
        if a >= e:
            break
        if b > s:
            yield lst[i]
        i += 1


rows = []
with open(catalog_in) as fh:
    for r in csv.DictReader(fh):
        if r["layer"] != "ERV":
            continue
        chrom, s, e = r["chr"], int(r["start"]), int(r["end"])

        best = None
        for a, b, strand, gid, im, cat in hits(prov.get(chrom, []), s, e):
            ov = min(b, e) - max(a, s)
            if best is None or ov > best[0]:
                best = (ov, a, b, gid, im, cat)
        if best is None:
            # Fallback: use the rmsk fragment assignment directly
            sub = collections.Counter()
            for a, b, name, reg, _st in hits(frag.get(chrom, []), s, e):
                if reg == "internal":
                    sub[name] += min(b, e) - max(a, s)
            top = sub.most_common(1)[0][0] if sub else ""
            rows.append([r["ervmap_id"], chrom, s, e, r["strand"], "", "", top, "", "", "", 0, top, "rmsk_fallback"])
            continue

        ov, a, b, gid, im, cat = best
        # rmsk fragments within the locus: split into internal regions and LTRs, then assign 5'/3' by relative position
        ltrs = []
        for fa, fb, name, reg, st in hits(frag.get(chrom, []), a, b):
            if reg == "ltr":
                ltrs.append((fa, fb, name))
        ltrs.sort()
        ltr5 = ltrs[0][2] if ltrs else ""
        ltr3 = ltrs[-1][2] if len(ltrs) > 1 else ""
        rows.append([
            r["ervmap_id"], chrom, s, e, r["strand"],
            gid, cat, im, ltr5, ltr3,
            ";".join(sorted({x[2] for x in ltrs})), len(ltrs),
            im, "provirus",
        ])

with open(out, "w") as fo:
    fo.write("ervmap_id\tchr\tstart\tend\tstrand\tprovirus_locus\tprovirus_category\t"
             "int_model\tltr5_model\tltr3_model\tltr_models\tn_ltr\tfallback_subfamily\tsource\n")
    for row in rows:
        fo.write("\t".join(map(str, row)) + "\n")

print(f"rows={len(rows)}", file=sys.stderr)
