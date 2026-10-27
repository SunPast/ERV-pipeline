# Output format

## 1. Unified long table

`results/aggregate/erv_expression_long.tsv.gz`

| Column | Type | Description |
|---|---|---|
| tool | text | `ervmap` / `telescope` / `squire` / `tetranscripts` |
| cohort | text | cohort (consistent with `Accession` in the database, e.g. `EGAD00001003977`) |
| sample | text | platform-anonymized sample ID (`<Accession>_S####` / `_P####`) |
| level | text | `locus` (locus level) or `subfamily` (subfamily level) |
| feature | text | locus_id when level=locus; subfamily name when level=subfamily |
| region_type | text | `ltr` / `internal`; empty at locus level for SQuIRE and ERVmap |
| count | numeric | raw count (not normalized) |

Notes:

- Each (tool, cohort, sample, level, feature) should appear only once;
- locus_id from different tools **cannot be compared across tools** (Telescope/TEtranscripts use the same set; SQuIRE uses its own instances; ERVmap uses official locus numbers);
- Subfamily level can be compared across tools, but pay attention to each tool's coverage (550 / 545 / 145 / 60).

## 2. Subfamily evidence table

`results/aggregate/erv_subfamily_evidence.tsv.gz`

| Column | Description |
|---|---|
| cohort | cohort |
| subfamily | subfamily name (consistent with rmsk `repName`) |
| region_type | `ltr` / `internal` (determined by name: names ending in `-int` are internal) |
| n_tools_support | how many of TEtranscripts / SQuIRE / Telescope detect it in this cohort (count>0) |
| ervmap_support | whether ERVmap also has a value (0/1), for external validation only |

## 3. Annotation validation

`results/annotation/annotation_check.tsv`

| metric | Meaning |
|---|---|
| canonical_loci | number of fragments in the unified reference |
| canonical_subfamilies | number of subfamilies in the unified reference |
| squire_loci | number of fragments in the SQuIRE BED |
| subfamily_overlap_med | median per-subfamily locus overlap rate |
| subfamily_ge95 | number of subfamilies with overlap rate ≥95% |
| provirus_loci | number of provirus annotation loci |
| provirus_with_ltr_category | number of loci of type prototype/oneside (containing LTRs) |
| provirus_internal_models | number of internal-region models covered by the provirus annotation |

## 4. ERVmap provirus composition mapping

`results/annotation/ervmap_provirus_map.tsv`

| Column | Description |
|---|---|
| ervmap_id | ERVmap locus number (corresponds to column 4 of `ERVmap.bed`) |
| provirus_locus | matched provirus locus (empty if no match) |
| provirus_category | prototype / oneside / internal |
| int_model | internal-region model (primary assignment) |
| ltr5_model / ltr3_model | rmsk names of the two flanking LTRs |
| ltr_models | all LTR names within the locus (semicolon-separated) |
| n_ltr | number of LTR elements |
| fallback_subfamily | rmsk assignment when no provirus matches |
| source | `provirus` / `rmsk_fallback` |
