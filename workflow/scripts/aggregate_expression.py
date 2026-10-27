#!/usr/bin/env python3
"""Aggregate the tools' outputs into a unified long table + subfamily evidence table.

Usage: aggregate_expression.py OUTDIR SAMPLES_TSV LOCUS_MAP ERVMAP_MAP MIN_COUNT

Unified long table columns:
  tool, cohort, sample, level(locus|subfamily), feature, region_type, count, milli_div

milli_div (RepeatMasker divergence, per mille) is a property of the locus, not of
the tool, so it is filled at locus level only, by joining the locus map. It is
left empty at subfamily level: a subfamily spans a wide divergence range, so a
single summary number would be misleading. ERVmap loci are curated proviruses,
not rmsk copies, so they get no milli_div either.

Subfamily evidence table columns:
  cohort, subfamily, region_type, n_tools_support, ervmap_support

Parsing notes (all verified against real output files):
  * Telescope report : col0 = locus id, col2 = final_count; first line is "## RunInfo".
  * TEcount cntTable : "<locus>:<Family>:<Class>\\t<count>"; gene rows are quoted ids.
  * SQuIRE TEcounts  : col3 = TE_ID (chr|start|end|Subfamily:Family:Class|milliDiv|strand),
                       col15 = tot_counts.
  * SQuIRE subFcounts: col2 = Subfamily:Family:Class, col6 = tot_counts.
  * ERVmap           : normalized by normalize_ervmap.py to locus_id / count.

Only ERV features (those present in the locus map) are written; SQuIRE also
reports non-ERV repeat subfamilies and those are dropped here on purpose.
"""
import collections
import csv
import gzip
import os
import sys

outdir, samples_tsv, locus_map, ervmap_map, min_count = sys.argv[1:6]
min_count = float(min_count)

samples = list(csv.DictReader(open(samples_tsv), delimiter="\t"))

locus2sf = {}
locus_md = {}
erv_subfamilies = set()
for r in csv.DictReader(open(locus_map), delimiter="\t"):
    locus2sf[r["locus_id"]] = (r["subfamily"], r["region_type"])
    locus_md[r["locus_id"]] = r.get("milli_div", "")
    erv_subfamilies.add(r["subfamily"])

erv2sf = {}
for r in csv.DictReader(open(ervmap_map), delimiter="\t"):
    sf = r["int_model"] or r["fallback_subfamily"]
    erv2sf[r["ervmap_id"]] = sf

long_rows = []


def add(tool, cohort, sample, level, feature, reg, count, md=""):
    long_rows.append((tool, cohort, sample, level, feature, reg, count, md))


def squire_locus_id(te_id):
    """chr|start|end|Subfamily:Family:Class|milliDiv|strand -> canonical locus_id"""
    p = te_id.split("|")
    if len(p) < 6:
        return None
    chrom, start, end, subf, _md, strand = p[0], p[1], p[2], p[3], p[4], p[5]
    return f"{subf.split(':')[0]}__{chrom}_{start}_{end}_{strand}"


for s in samples:
    cohort, sample = s["cohort"], s["sample"]
    # subfamily totals are summed over loci; emitting one row per locus here
    # would duplicate the same subfamily thousands of times
    subsum = collections.defaultdict(float)

    # ---- Telescope (locus) ----
    rep = f"{outdir}/telescope/{cohort}/{sample}/{sample}.report.tsv"
    if os.path.exists(rep):
        with open(rep) as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("transcript\t"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                gid, cnt = f[0], float(f[2])
                if gid not in locus2sf:
                    continue
                sf, reg = locus2sf[gid]
                add("telescope", cohort, sample, "locus", gid, reg, cnt, locus_md.get(gid, ""))
                subsum[("telescope", sf, reg)] += cnt

    # ---- TEtranscripts / TEcount (locus) ----
    cnt_tbl = f"{outdir}/tetranscripts/{cohort}/{sample}/{sample}.cntTable"
    if os.path.exists(cnt_tbl):
        with open(cnt_tbl) as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2:
                    continue
                name = f[0].strip('"').split(":")[0]
                if name not in locus2sf:
                    continue
                sf, reg = locus2sf[name]
                cnt = float(f[1])
                if cnt <= 0:
                    continue
                add("tetranscripts", cohort, sample, "locus", name, reg, cnt, locus_md.get(name, ""))
                subsum[("tetranscripts", sf, reg)] += cnt

    # ---- SQuIRE locus level ----
    sq_te = f"{outdir}/squire/{cohort}/{sample}/{sample}_TEcounts.txt"
    if os.path.exists(sq_te):
        with open(sq_te) as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 16:
                    continue
                lid = squire_locus_id(f[3])
                if lid not in locus2sf:
                    continue
                sf, reg = locus2sf[lid]
                add("squire", cohort, sample, "locus", lid, reg, float(f[15]), locus_md.get(lid, ""))

    # ---- SQuIRE subfamily level ----
    sq_sf = f"{outdir}/squire/{cohort}/{sample}/{sample}_subFcounts.txt"
    if os.path.exists(sq_sf):
        with open(sq_sf) as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 7:
                    continue
                sf = f[2].split(":")[0]
                if sf not in erv_subfamilies:
                    continue
                reg = "internal" if sf.endswith("-int") else "ltr"
                add("squire", cohort, sample, "subfamily", sf, reg, float(f[6]))

    # ---- ERVmap (locus; subfamily from the mapping table) ----
    ev = f"{outdir}/ervmap/{cohort}/{sample}/{sample}.locus.tsv"
    if os.path.exists(ev):
        with open(ev) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                cnt = float(r["count"])
                add("ervmap", cohort, sample, "locus", r["locus_id"], "", cnt)
                sf = erv2sf.get(r["locus_id"], r["subfamily"])
                if sf:
                    reg = "internal" if sf.endswith("-int") else "ltr"
                    subsum[("ervmap", sf, reg)] += cnt

    for (tool, sf, reg), v in subsum.items():
        add(tool, cohort, sample, "subfamily", sf, reg, v)

# ---- Write long table ----
os.makedirs(f"{outdir}/aggregate", exist_ok=True)
with gzip.open(f"{outdir}/aggregate/erv_expression_long.tsv.gz", "wt") as fo:
    fo.write("tool\tcohort\tsample\tlevel\tfeature\tregion_type\tcount\tmilli_div\n")
    for row in long_rows:
        fo.write("\t".join(map(str, row)) + "\n")

# ---- Per-sample QC: a tool that produced nothing must be visible, not silent ----
qc = collections.defaultdict(lambda: [0.0, 0])
for tool, cohort, sample, level, feature, reg, cnt, _md in long_rows:
    if level != "locus":
        continue
    k = (tool, cohort, sample)
    qc[k][0] += cnt
    if cnt > 0:
        qc[k][1] += 1
with open(f"{outdir}/aggregate/qc_per_sample.tsv", "w") as fq:
    fq.write("tool\tcohort\tsample\ttotal_counts\tn_positive_loci\tflag\n")
    for (tool, cohort, sample), (tot, npos) in sorted(qc.items()):
        flag = "ZERO_COUNTS" if tot <= 0 else ("LOW_COUNTS" if tot < 1000 else "OK")
        fq.write(f"{tool}\t{cohort}\t{sample}\t{tot:.1f}\t{npos}\t{flag}\n")

# ---- Subfamily evidence table ----
# A tool supports a subfamily in a cohort if any sample reaches min_count.
sup = collections.defaultdict(set)
for tool, cohort, sample, level, feature, reg, cnt, _md in long_rows:
    if level != "subfamily" or cnt < min_count:
        continue
    sup[(cohort, feature)].add(tool)

agg = collections.defaultdict(lambda: [0, 0, ""])
for (cohort, sf), tools in sup.items():
    n = len(tools - {"ervmap"})
    e = 1 if "ervmap" in tools else 0
    reg = "internal" if sf.endswith("-int") else "ltr"
    agg[(cohort, sf)] = [n, e, reg]

with gzip.open(f"{outdir}/aggregate/erv_subfamily_evidence.tsv.gz", "wt") as fo:
    fo.write("cohort\tsubfamily\tregion_type\tn_tools_support\tervmap_support\n")
    for (cohort, sf), (n, e, reg) in sorted(agg.items()):
        fo.write(f"{cohort}\t{sf}\t{reg}\t{n}\t{e}\n")

print(f"long_rows={len(long_rows)} evidence={len(agg)}", file=sys.stderr)
