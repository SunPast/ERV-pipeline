#!/usr/bin/env python3
"""Export the per-tool "standard" CSVs consumed by the database loader.

The SQL side (ERV-data/rerun_out/db_load/build_import_csvs.sh) does not read the
unified long table; it reads one CSV per tool/level. This step produces exactly
those files, so the loading scripts need no changes apart from the file names.

Usage:
  export_standard.py OUTDIR SAMPLES_TSV LOCUS_MAP ERVMAP_MAP STANDARD_DIR

Outputs (all comma separated, header on line 1, zero-count rows removed):

  ervmap_standard.csv
      cohort,sample,chr,start,end,ervmap_id,strand,count
  tetranscripts_subf_standard.csv
      cohort,sample,subfamily,family,count
  tetranscripts_locus_standard.csv
      cohort,sample,locus,subfamily,family,chr,start,end,strand,count
  telescope_locus_standard.csv
      cohort,sample,locus,subfamily,family,chr,start,end,strand,count
  squire_te_standard.csv
      cohort,sample,subfamily,family,count,TE_chr,TE_start,TE_stop,milliDiv,TE_strand
  squire_subf_standard.csv
      cohort,sample,subfamily,family,count

Notes:
  * sample names stay raw; the SQL side maps them to pseudonyms and adds
    norm_count from its own size factors.
  * only ERV loci/subfamilies present in the locus map are written; SQuIRE's
    non-ERV reference rows are dropped here.
  * the loader keeps only families ERV1/ERVK/ERVL/ERVL-MaLR, so the 3,387 loci
    with family ERVL? / ERV1? are accepted to be dropped at load time.
"""
import csv
import os
import sys

outdir, samples_tsv, locus_map, ervmap_map, stddir = sys.argv[1:6]

samples = list(csv.DictReader(open(samples_tsv), delimiter="\t"))

meta = {}            # locus_id -> (subfamily, family, chr, start, end, strand)
erv_subfamilies = set()
for r in csv.DictReader(open(locus_map), delimiter="\t"):
    meta[r["locus_id"]] = (
        r["subfamily"], r["family"], r["chr"], r["start"], r["end"], r["strand"],
        r.get("milli_div", ""),
    )
    erv_subfamilies.add(r["subfamily"])

erv2sf = {}
for r in csv.DictReader(open(ervmap_map), delimiter="\t"):
    erv2sf[r["ervmap_id"]] = r["int_model"] or r["fallback_subfamily"]


def squire_locus_id(te_id):
    p = te_id.split("|")
    if len(p) < 6:
        return None
    return f"{p[3].split(':')[0]}__{p[0]}_{p[1]}_{p[2]}_{p[5]}"


os.makedirs(stddir, exist_ok=True)
paths = {
    "ervmap": f"{stddir}/ervmap_standard.csv",
    "tt_subf": f"{stddir}/tetranscripts_subf_standard.csv",
    "tt_locus": f"{stddir}/tetranscripts_locus_standard.csv",
    "tel_locus": f"{stddir}/telescope_locus_standard.csv",
    "sq_te": f"{stddir}/squire_te_standard.csv",
    "sq_subf": f"{stddir}/squire_subf_standard.csv",
}
fh = {k: open(v, "w", newline="") for k, v in paths.items()}
w = {k: csv.writer(fh[k]) for k in fh}

w["ervmap"].writerow(["cohort", "sample", "chr", "start", "end", "ervmap_id", "strand", "count"])
w["tt_subf"].writerow(["cohort", "sample", "subfamily", "family", "count"])
w["tt_locus"].writerow(["cohort", "sample", "locus", "subfamily", "family", "chr", "start", "end", "strand", "count", "milli_div"])
w["tel_locus"].writerow(["cohort", "sample", "locus", "subfamily", "family", "chr", "start", "end", "strand", "count", "milli_div"])
w["sq_te"].writerow(["cohort", "sample", "subfamily", "family", "count", "TE_chr", "TE_start", "TE_stop", "milliDiv", "TE_strand"])
w["sq_subf"].writerow(["cohort", "sample", "subfamily", "family", "count"])

for s in samples:
    cohort, sample = s["cohort"], s["sample"]

    # ---- ERVmap ----
    ev = f"{outdir}/ervmap/{cohort}/{sample}/{sample}.locus.tsv"
    if os.path.exists(ev):
        for r in csv.DictReader(open(ev), delimiter="\t"):
            cnt = float(r["count"])
            if cnt <= 0:
                continue
            w["ervmap"].writerow([cohort, sample, r["chr"], r["start"], r["end"],
                                  r["locus_id"], r["strand"], cnt])

    # ---- TEtranscripts (locus + subfamily) ----
    tt = f"{outdir}/tetranscripts/{cohort}/{sample}/{sample}.cntTable"
    subf = {}
    if os.path.exists(tt):
        with open(tt) as f:
            next(f, None)
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 2:
                    continue
                locus = p[0].strip('"').split(":")[0]
                if locus not in meta:
                    continue
                cnt = float(p[1])
                if cnt <= 0:
                    continue
                sf, fam, chrom, start, end, strand, md = meta[locus]
                w["tt_locus"].writerow([cohort, sample, locus, sf, fam,
                                        chrom, start, end, strand, cnt, md])
                # subfamily counts are sums over loci, not one row per locus
                key = (sf, fam)
                subf[key] = subf.get(key, 0.0) + cnt

        for (sf, fam), v in subf.items():
            w["tt_subf"].writerow([cohort, sample, sf, fam, v])

    # ---- Telescope ----
    tel = f"{outdir}/telescope/{cohort}/{sample}/{sample}.report.tsv"
    if os.path.exists(tel):
        with open(tel) as f:
            for line in f:
                if line.startswith("#") or line.startswith("transcript\t"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 3:
                    continue
                locus = p[0]
                if locus not in meta:
                    continue
                cnt = float(p[2])
                if cnt <= 0:
                    continue
                sf, fam, chrom, start, end, strand, md = meta[locus]
                w["tel_locus"].writerow([cohort, sample, locus, sf, fam,
                                         chrom, start, end, strand, cnt, md])

    # ---- SQuIRE locus ----
    sq = f"{outdir}/squire/{cohort}/{sample}/{sample}_TEcounts.txt"
    if os.path.exists(sq):
        with open(sq) as f:
            next(f, None)
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 16:
                    continue
                lid = squire_locus_id(p[3])
                if lid not in meta:
                    continue
                cnt = float(p[15])
                if cnt <= 0:
                    continue
                sf, fam = meta[lid][0], meta[lid][1]
                w["sq_te"].writerow([cohort, sample, sf, fam, cnt,
                                     p[8], p[9], p[10], p[12], p[13]])

    # ---- SQuIRE subfamily ----
    sqs = f"{outdir}/squire/{cohort}/{sample}/{sample}_subFcounts.txt"
    if os.path.exists(sqs):
        with open(sqs) as f:
            next(f, None)
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 7:
                    continue
                sf, fam = p[2].split(":")[0], p[2].split(":")[1]
                if sf not in erv_subfamilies:
                    continue
                cnt = float(p[6])
                if cnt <= 0:
                    continue
                w["sq_subf"].writerow([cohort, sample, sf, fam, cnt])

for f in fh.values():
    f.close()

print("standard tables written to " + stddir, file=sys.stderr)
