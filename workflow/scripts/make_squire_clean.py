#!/usr/bin/env python3
"""Generate SQuIRE's clean_folder from te_locus_map.tsv (two files, no Build/Clean run).

Outputs (column formats copied exactly from native SQuIRE squeue_clean output, see $squire_ref/squire_clean/):
  hg38_all.bed         9 columns: chrom start end TE_ID milliDiv strand start end rgb
                       TE_ID = chrom|start|end|Subfamily:Family:Class|milliDiv|strand
  hg38_all_copies.txt  4 columns: Subfamily:Family:Class  Copies  Tot.Length  Avg.Length

Why build it ourselves:
  squire Build automatically downloads the current UCSC rmsk over the network and uses bedtools merge to merge adjacent fragments.
  What we want is "not merged, same coordinates as TEtranscripts/Telescope", so we bypass Build and
  lay out the unified annotation in a form Count can read (Count only needs one *.bed + one *_copies.txt in the directory).

milliDiv: RepeatMasker's divergence (per mille from the consensus model).
          The TEtranscripts GTF does not carry it, so it is looked up by
          (chrom, start, end, repName) in the UCSC rmsk table when that table is
          supplied as the third argument; otherwise it falls back to 0.
          It does not affect counting (Count.py only writes it out; only the
          unused `squire Call` module reads it), but it is kept so the SQuIRE
          output columns stay meaningful.

Only primary chromosomes (chr1-22, X, Y) are written, for the same reason as in
build_locus_gtf.py: alt / fix / random / chrUn names never match the contig names
in this project's BAM headers (bare accessions such as GL000008.2).
"""
import csv
import re
import sys
from collections import defaultdict

map_in, outdir = sys.argv[1:3]
rmsk_table = sys.argv[3] if len(sys.argv) > 3 else None

PRIMARY = re.compile(r"^chr([0-9]+|X|Y)$")

millidiv = {}
if rmsk_table:
    with open(rmsk_table) as fh:
        for line in fh:
            p = line.split()
            if len(p) < 17:
                continue
            millidiv[(p[5], p[6], p[7], p[10])] = p[2]

copies = defaultdict(int)
lengths = defaultdict(int)

with open(f"{outdir}/hg38_all.bed", "w") as fb:
    with open(map_in) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            chrom = r["chr"]
            if not PRIMARY.match(chrom):
                continue
            start = int(r["start"])
            end = int(r["end"])
            strand = r["strand"]
            sub = f'{r["subfamily"]}:{r["family"]}:{r["rep_class"]}'
            md = r.get("milli_div") or millidiv.get((chrom, str(start), str(end), r["subfamily"]), "0")
            te_id = f"{chrom}|{start}|{end}|{sub}|{md}|{strand}"
            fb.write(
                "\t".join(
                    [chrom, str(start), str(end), te_id, md, strand,
                     str(start), str(end), "94,137,255"]
                ) + "\n"
            )
            copies[sub] += 1
            lengths[sub] += end - start

with open(f"{outdir}/hg38_all_copies.txt", "w") as fc:
    fc.write("Subfamily:Family:Class\tCopies\tTot.Length\tAvg.Length\n")
    for sub in sorted(copies):
        fc.write(f"{sub}\t{copies[sub]}\t{lengths[sub]}\t{lengths[sub]/copies[sub]:.2f}\n")

print(f"bed_rows={sum(copies.values())} subfamilies={len(copies)}", file=sys.stderr)
