# Output format

All paths below are relative to `outdir` in `config/config.yaml`
(node-local disk, e.g. `/data3/dingjia/erv_rerun_out/pipeline_results`).

## 1. Unified long table

`aggregate/erv_expression_long.tsv.gz`

| Column | Type | Description |
|---|---|---|
| tool | text | `ervmap` / `telescope` / `squire` / `tetranscripts` |
| cohort | text | cohort label (consistent with the database, e.g. `EGAD00001003977`) |
| sample | text | sample ID as in the manifest (pseudonymisation happens at load time) |
| level | text | `locus` or `subfamily` |
| feature | text | locus_id when level=locus; subfamily name when level=subfamily |
| region_type | text | `ltr` / `internal` |
| count | numeric | raw count (not normalized; FPKM is never used) |
| milli_div | integer | RepeatMasker divergence, per mille. Filled at **locus level only** (never for ERVmap); empty at subfamily level |

Notes:

- Each (tool, cohort, sample, level, feature) should appear only once.
- Telescope, TEtranscripts and SQuIRE share the **same 693,078 locus ids**
  (byte-identical annotation), so locus level is directly comparable between those three.
  ERVmap keeps its own 3,220 provirus ids and only joins at the subfamily level.
- Subfamily level is comparable across all four tools; coverage is 550 subfamilies for the
  three rmsk-based tools and 145 for ERVmap.

## 2. Subfamily evidence table

`aggregate/erv_subfamily_evidence.tsv.gz`

| Column | Description |
|---|---|
| cohort | cohort |
| subfamily | subfamily name (consistent with rmsk `repName`) |
| region_type | `ltr` / `internal` (names ending in `-int` are internal) |
| n_tools_support | how many of TEtranscripts / SQuIRE / Telescope reach `min_count` in this cohort |
| ervmap_support | 1/0, whether ERVmap also supports it (external validation only) |

## 3. Annotation validation

`annotation/annotation_check.tsv`

| metric | Meaning |
|---|---|
| canonical_loci | number of fragments in the unified reference |
| canonical_subfamilies | number of subfamilies in the unified reference |
| squire_loci | number of fragments in the SQuIRE BED |
| subfamily_overlap_med | median per-subfamily locus overlap rate |
| subfamily_ge95 | number of subfamilies with overlap rate >= 95% |
| provirus_loci | number of provirus annotation loci |
| provirus_with_ltr_category | number of prototype/oneside loci (containing LTRs) |
| provirus_internal_models | number of internal-region models covered |

## 4. ERVmap provirus composition mapping

`annotation/ervmap_provirus_map.tsv`

| Column | Description |
|---|---|
| ervmap_id | ERVmap locus number (column 4 of `ERVmap.bed`) |
| provirus_locus | matched provirus locus (empty if no match) |
| provirus_category | prototype / oneside / internal |
| int_model | internal-region model (primary assignment) |
| ltr5_model / ltr3_model | rmsk names of the two flanking LTRs |
| ltr_models | all LTR names within the locus (semicolon separated) |
| n_ltr | number of LTR elements |
| fallback_subfamily | rmsk assignment when no provirus matches |
| source | `provirus` / `rmsk_fallback` |

## 5. Standard tables for the database loader

`standard/` — one CSV per tool/level, comma separated, header on line 1, zero-count rows
removed, sample names left raw (the loader applies pseudonyms and `norm_count`).

| File | Columns | Target table |
|---|---|---|
| `ervmap_standard.csv` | cohort, sample, chr, start, end, ervmap_id, strand, count | tool_ervmap |
| `tetranscripts_subf_standard.csv` | cohort, sample, subfamily, family, count | tool_tetranscripts |
| `tetranscripts_locus_standard.csv` | cohort, sample, locus, subfamily, family, chr, start, end, strand, count, milli_div | tool_tetranscripts_locus |
| `telescope_locus_standard.csv` | cohort, sample, locus, subfamily, family, chr, start, end, strand, count, milli_div | tool_telescope |
| `squire_te_standard.csv` | cohort, sample, subfamily, family, count, TE_chr, TE_start, TE_stop, milliDiv, TE_strand | tool_squire_te |
| `squire_subf_standard.csv` | cohort, sample, subfamily, family, count | tool_squire_subf |
