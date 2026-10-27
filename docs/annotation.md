# Annotation: sources, differences, and decisions

This document answers a recurring question: **do the annotations of the four tools actually need to be unified?**

## 1. What the four annotations are

| Tool | Annotation | Source | Form | Coverage |
|---|---|---|---|---|
| TEtranscripts | `hg38_rmsk_TE.gtf` | Released by the Hammell lab with TEtranscripts | rmsk fragments, `gene_id`=repName | 550 subfamilies / 720,871 fragments |
| Telescope (old) | `telescope_full_plus23.gtf` | `mlbendall/telescope_annotation_db` → `HERV_rmsk.hg38.v2` (RepeatMasker) | provirus loci (fragments merged), 14,968 entries | 60 internal-region models |
| SQuIRE | `hg38_all.bed` | Downloaded automatically by local `squire Build` from UCSC rmsk (2025-04) | merged TE instances | 545 subfamilies / 702,850 instances |
| ERVmap | `ERVmap.bed` | Official repository `ref/ERVmap.bed` from Tokuyama et al. 2018 | curated provirus loci | 3,220 loci / 145 subfamilies |

Additional facts:

- Telescope's official annotation **does not contain LTR element names**, but the locus span **does include LTRs** (prototype 100% contain 2 LTRs, oneside 100% contain 1).
- Column 6 (strand) of `ERVmap.bed` is **all `+`**; the real strand orientation is in `final/locus_catalog.csv` (1,625 `+` / 1,587 `-` / 8 `.`).
- The locus overlap between SQuIRE and TEtranscripts has a per-subfamily median of **96.8%**; 545/550 are above 80%, the only outlier being MER70B (45%).

## 2. Why we do not force the inputs to be unified

- **ERVmap**: the BED is the method itself. With rmsk fragments, the premise of "strict best-match + mismatch filtering" no longer holds, and the result is no longer ERVmap.
- **SQuIRE 0.9.9.92**: only `squire Build` (auto-downloads rmsk) is available, there is no `MakeBed`, and it does not accept external GTF/BED.
- **Telescope**: accepts any valid GTF, so it **can** use the unified reference; it is also the only tool that **must be rerun** (the old annotation had only 60 models).
- **TEtranscripts**: it is itself built on the rmsk fragment system, so using it as the unified reference is the most natural choice.

## 3. Decisions

| Decision | Content |
|---|---|
| Reference | Official TEtranscripts rmsk GTF (550 subfamilies) |
| Derived | One locus GTF (720,871, fragment level, not merged) shared by TEtranscripts and Telescope |
| SQuIRE | Keep native annotation, no rerun, no mapping; aligned at subfamily level, locus level not compared |
| ERVmap | No rerun; only map "provirus locus → internal region/LTR composition" |

## 4. Why not use "2 out of 4"

The four tools cover different spaces (60 / 145 / 545 / 550); a 2/4 rule would effectively let ERVmap decide inclusion.
Instead we use a **layered consensus**:

- LTR layer: TEtranscripts + SQuIRE + Telescope, 2/3
- Internal-region layer: same as above, 2/3
- Provirus layer: ERVmap + Telescope (old 60), 2/2, a separate layer
- Four-tool common layer: 60 internal regions, display only, no threshold

ERVmap's role is "external independent validation", recorded separately in the `ervmap_support` field.

## 5. Provirus composition mapping rules (ERVmap)

1. For each ERVmap locus, find the provirus locus with the largest overlap;
2. Primary assignment = that locus's `intModel` (internal region);
3. LTR names are taken directly from the rmsk fragments inside the locus (`ltr5_model` / `ltr3_model`);
4. **No coverage-fraction threshold**: the internal-region coverage distribution is continuous with no valley, so any threshold would be arbitrary;
5. Loci that do not match a provirus locus (measured: 399) fall back to the rmsk fragment assignment.

**LTR expression is not derived from ERVmap**: a provirus is one transcript, counting each of its two LTRs would double-count, and solo LTRs would be missed.
LTR expression uses the measured values from TEtranscripts / SQuIRE.

## 6. Known biases and how they are handled

| Bias | Value | Handling |
|---|---|---|
| SQuIRE vs reference locus difference | 3.2% | Comparable at subfamily level, not compared at locus level; noted in the report |
| Old Telescope annotation only reached internal regions | 60 vs 550 | Rerun with the unified GTF |
| ERVmap has only 145 subfamilies | 87 internal + 58 LTR | External validation only, does not participate in the consensus threshold |

## Single-version source: the UCSC rmsk table

Because milliDiv (RepeatMasker divergence) is used by downstream analysis,
coordinates and milliDiv must come from the **same** library release. The
TEtranscripts GTF has no milliDiv column, so the annotation is now built
directly from the UCSC rmsk table that SQuIRE downloaded (2025-04):

```
input : /data/data2/reference/hg38_squire2/hg38_rmsk.txt   (5,683,690 rows)
script: workflow/scripts/build_locus_from_rmsk.py
filter: chr1-22,X,Y  and  the 550 ERV subfamilies from the dictionary
output: 693,078 loci, each with subfamily / family / class / strand / milliDiv
```

Verification that this changes nothing about the loci:

```
rmsk-derived loci      : 693,078
GTF-derived loci       : 693,078
differences on (chr,start,end,strand,subfamily): 0
```

So the switch removes the cross-version milliDiv problem at zero cost: the
locus set is byte-identical to the one previously derived from the
TEtranscripts GTF. Only 192 of 693,078 loci report milliDiv = 0.

Consequences:
- `te_locus_map.tsv` gains a `milli_div` column (metadata for the frontend).
- the SQuIRE bed takes milliDiv straight from the map (same release).
- the TEtranscripts GTF is no longer needed for coordinates; it stays only as
  the historical provenance of the ERV subfamily dictionary.

## Chromosome scope: primary assembly only

Decision (2026-10-27): keep `chr1-22, X, Y` only.

| Scope | ERV loci |
|---|---|
| all contigs | 720,871 |
| primary only | **693,078** |
| dropped (alt / fix / random / chrUn) | 27,793 |

Reasons, all measured:

1. **Contig names do not match the BAMs.** BAM headers contain 194 sequences:
   chr1-22/X/Y/M plus 169 bare accessions (`GL000008.2`, `KI270706.1`, ...).
   They contain **zero** sequences named `*_alt`, `*_fix`, `*_random` or `chrUn*`,
   so every annotation row on those contigs receives 0 reads. Of the 27,793
   dropped loci only 2,436 can even be renamed to a contig present in the BAM.
2. **Reads there are mostly multi-mappers.** On the three largest matched contigs,
   52-90% of reads have NH>1: the same molecule also maps to the primary copy.
3. **No subfamily-level loss.** All 259 subfamilies involved in those 2,436 loci
   also have loci on primary chromosomes.
4. **chrM carries no ERV.** rmsk has 5 rows on chrM (rRNA / tRNA / simple repeat)
   and neither annotation includes them. chrM reads (0.59% in the test sample)
   must be excluded from library-size / size-factor estimation.

Verification after the filter:

```
our ERV loci on primary   : 693078
SQuIRE ERV loci on primary: 693078
differences               : 0
```

The two ERV annotations are now identical on the primary assembly. The only
remaining difference is that SQuIRE's file also carries non-ERV repeat classes,
which it needs as a multi-mapper reference and which are dropped at aggregation.

## SQuIRE measurements (2026-10-27, four Count runs on the same BAM, only clean_folder changed)

Sample: `GO29293_ngs_rna_targrna_SAM04c589eb3fb3_20160128` (IMvigor210), first 1M read pairs,
`--read_length 36 --pthreads 10`, measured aligned_libsize = 978,247.
**All comparisons below are restricted to the 545 shared ERV subfamilies** (positive counts from different universes are not comparable).

| clean_folder | Subfamilies | Total ERV counts | Positive ≥1 / ≥5 / ≥10 | MSTA | r vs native |
|---|---|---|---|---|---|
| Native (all repeat classes, rmsk 2025-04) | 1,280 | 1,961 | 221 / 65 / 33 | 622.22 | — |
| ERV only (550 / 720,871 fragments) | 550 | 1,305 | 220 / 58 / 28 | 57.00 | 0.5469 |
| Union (our ERV + native non-ERV) | 1,286 | 1,961 | 221 / 65 / 33 | 622.22 | 1.00000 |
| **All TE (official TEtranscripts GTF)** | **1,180** | **2,011** | **221 / 65 / 34** | **623.08** | **0.99750** |

- Looking only at "detected or not", the four annotations are almost indistinguishable (221 vs 220); the difference is in **total amount**, concentrated in high-copy subfamilies
  (MSTA 622 → 57; its `uniq_counts` is 57 on both sides, and what is lost is the part EM assigned to multi-mapping reads)
- Hence the rule: **SQuIRE's annotation must be "all TE", but the ERV part can use our coordinates**

## The official TEtranscripts GTF is actually all-TE

`TEtranscripts_hg38_rmsk_TE.gtf.gz`: 4,693,511 loci / 1,180 subfamilies

| class | Entries |
|---|---|
| SINE | 1,844,056 |
| LINE | 1,570,523 |
| LTR | 748,597 |
| DNA | 499,108 |
| Other (Satellite/Retroposon/Unknown/RC/RNA) | 31,227 |

The earlier `build_locus_gtf.py` filtered by the 550 ERV subfamilies in the dictionary, leaving only 720,871
ERV fragments — **that was caused by our filtering, the file itself is all-TE**.

The SQuIRE clean_folder generated from the full GTF is in
`/data/home/dingjia/ERV-data/annotation/squire_clean_full/` (4,693,511-row bed + copies for 1,180 subfamilies).
