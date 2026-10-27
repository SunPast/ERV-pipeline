# Work order: full ERV/TE quantification rerun

Status: ready to run, pending the 4,289-sample manifest.
Owner: backend. Frontend is not required for this work order.

## 0. Node and storage policy (mandatory)

- Execute on **node3 only**. No step may run on the master node.
- Working copy: `/data3/dingjia/ERV-pipeline_run` (rsync of the archived code).
- Products: `/data3/dingjia/erv_rerun_out/pipeline_results` (node3 local disk, /dev/sda).
- The archived code lives in `/data/home/dingjia/ERV-pipeline` and is the
  reproducibility source. It contains: config, workflow rules, scripts, docs.
  Re-create a run with:
  ```
  rsync -a --exclude .venv --exclude tests /data/home/dingjia/ERV-pipeline/ /data3/dingjia/ERV-pipeline_run/
  cd /data3/dingjia/ERV-pipeline_run && /data/home/dingjia/ERV-pipeline/.venv/bin/snakemake -p --cores 160 all --config samples=config/samples.tsv
  ```
- Inputs are read-only from `/data/data2` (BAM) and are not copied.
- After the run is reviewed: migrate the three aggregate tables (and, if wanted,
  per-sample outputs) back to `/data/home/dingjia/ERV-data/rerun_out_new/`,
  then delete superseded old products.

## 1. What is already verified (do not re-argue)

| Item | Value | How it was verified |
|---|---|---|
| Annotation source | official TEtranscripts hg38 RepeatMasker GTF | file `TEtranscripts_hg38_rmsk_TE.gtf.gz` |
| ERV loci (primary chromosomes, no merging) | **693,078** / 550 subfamilies | `wc -l te_locus_map.tsv` minus header |
| Full-TE reference for SQuIRE | 4,520,928 loci / 1,180 subfamilies | `wc -l squire_clean_full/hg38_all.bed` |
| Our ERV loci vs SQuIRE native ERV loci, primary only | 693,078 vs 693,078, **0 differences** | row-wise `comm` on (chr, start, end, strand, subfamily) |
| SQuIRE needs full-TE reference, not ERV only | ERV-only input drops high-copy subfamilies (MSTA 622 -> 57) | 4 Count runs on one BAM, see docs/annotation.md |
| SQuIRE Count can read a plain STAR BAM (no Map) | exit 0, both outputs written, 3.2 GB, ~19 min/sample | test run |
| Telescope with 693k-locus GTF | assignment rate 96.7% (old 60-model run: 96.1%) | test run, 15.7M fragments |
| TEcount peak memory / time | 2.61 GB / 1:24:14 per sample | `/usr/bin/time -v`, test run |
| Alt / fix / random contigs | dropped: names never match BAM headers; include only 0.6% of reads that are mostly multimappers | see docs/annotation.md |

## 1. Phase 0 - sample manifest (blocking)

Target file: `config/samples.tsv`, columns `cohort, sample, run, bam, fq1, fq2`.

- 30 cohorts, 4,289 samples (source of truth: `/data/home/dingjia/ERV-data/rerun_out/analysis/sample_cohort_map.csv`).
- `bam` must point to the existing STAR BAM under `/data/data2/IO_RNA/AS/<cohort_dir>/aligned/<run>/`.
- `fq1/fq2` are only needed if `tools.squire.map: true`; they can be left empty otherwise.

Checkpoint: every row's `bam` exists; count = 4,289; no duplicate `sample`.

## 2. Phase 1 - annotation layer (runs once, ~2 min)

```
snakemake -p --cores 8 results/annotation/annotation_check.tsv
```

Outputs: `results/annotation/{te_locus.gtf, te_locus_map.tsv, te_locus.exon.gtf,
te_locus_map_full.tsv, squire_clean/*, annotation_check.tsv}`

Acceptance: ERV locus count = 693,078; full map = 4,520,928; SQuIRE bed = 4,520,928.

## 3. Phase 2 - Telescope (all samples)

- Per sample: `samtools collate` (16 threads) + `telescope assign --ncpu 1`.
- Output: `results/telescope/<cohort>/<sample>/<sample>.report.tsv`
- Measured: 16:42 for 15.7M fragments; projected 51 min for a median 47.7M sample.
- Projected total: ~6,900 core-hours -> **43 h at 160 cores**.
- Checkpoint after the first 20 samples: assignment rate (unique+ambig)/total >= 90%;
  report row count > 3,000; `< 50 bp` loci must carry < 10% of counts.

## 4. Phase 3 - TEtranscripts locus level (all samples)

- Command: `TEcount --sortByPos --format BAM --mode multi -b <bam> --GTF gencode.v34 --TE te_locus.exon.gtf`
- Output: `results/tetranscripts/<cohort>/<sample>/<sample>.cntTable`
- Measured: 1:24:14, 2.61 GB, 1 core. Projected **~38 h at 160 cores**.
- Checkpoint: cntTable has ~4.5M rows (all TE loci + genes); ERV rows are prefixed with locus ids.

## 5. Phase 4 - SQuIRE Count (all samples, no Map)

- Command: `squire Count --map_folder <out> --clean_folder results/annotation/squire_clean`
- Output: `<sample>_TEcounts.txt` (locus) + `<sample>_subFcounts.txt` (1,180 subfamilies; keep the 550 ERV ones)
- Measured: 1:36 for 1M read pairs, 3.2 GB. Projected ~19 min/sample -> **~11 h at 160 cores**.
- Hard requirement: the map folder must contain no `*.log`, otherwise `getlibsize()` crashes.
- Checkpoint: `_subFcounts.txt` has 1,180 rows; ERV subfamily totals are non-zero for >= 200 subfamilies.

## 6. Phase 5 - aggregation and QC

```
snakemake -p --cores 8 results/aggregate/erv_expression_long.tsv.gz
```

Outputs:
- `results/aggregate/erv_expression_long.tsv.gz` (tool, cohort, sample, level, feature, region_type, count)
- `results/aggregate/erv_subfamily_evidence.tsv.gz` (n_tools_support over the three rerun tools + ervmap_support)

Acceptance:
1. per-tool row counts are consistent with sample counts x feature counts;
2. subfamily-level Spearman between Telescope and TEtranscripts > 0.8 on a 50-sample subset;
3. no sample has zero counts for all three tools;
4. ERVmap subfamily support covers 145 subfamilies (its native scope).

## 7. Phase 6 - database load

- Load raw counts only. FPKM columns are not used anywhere.
- Normalisation: DESeq2 size factors computed on **genes only**; chrM reads excluded.
- `filter.min_count = 2` applies at subfamily level.

## 8. Resource summary at 160 cores

| Phase | Core-hours | Wall time |
|---|---|---|
| Telescope | ~6,900 | 43 h |
| TEcount | ~6,000 | 38 h |
| SQuIRE Count | ~1,400 | 11 h |
| Aggregation | < 5 | < 1 h |
| **Total** | **~14,300** | **~92 h (~3.8 days)** |

Memory is not a constraint after removing `squire Map` (peak per process: 3.2 GB).

## 9. Failure handling

- Every rule is restartable per sample; rerun the same target and finished samples are skipped.
- A sample that fails twice: record it in `results/qc_failed.tsv` and continue the cohort.
- Never mix old and new SQuIRE results in the database: the old tables must be replaced as a whole.

## 10. Open items

1. Generate the 4,289-row manifest (Phase 0).
2. Decide whether locus-level ERVmap import runs for all cohorts (currently driven by
   `results/ervmap_raw/<cohort>/<sample>ERVresults.txt`; missing files are skipped).
3. Write the normalisation note (genes only, chrM excluded) into the loader SQL before Phase 6.
