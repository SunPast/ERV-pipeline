# Parameter reference

## Global

| Parameter | Default | Description |
|---|---|---|
| `genome` | hg38 | Must match the BAM; do not mix with hg19 |
| `annotation_release` | GENCODE_v34 | Must match the annotation used when aligning the BAM |
| `filter.min_count` | 2 | Platform noise floor, consistent with the online database |
| `normalization.erv` | deseq2_poscounts | ERV-side normalization, reuses the existing size factors |

## TEtranscripts

| Parameter | Value | Description |
|---|---|---|
| `--mode` | multi | Multi-mapping reads are distributed by EM; `uniq` substantially underestimates |
| `--GTF` | GENCODE v34 | Gene annotation |
| `--TE` | unified locus GTF | One gene per rmsk fragment, which gives locus level |
| `level` | locus + subfamily | Subfamily is aggregated from locus, no separate run needed |
| Counting backend | `TEcount` | Single-sample counting; `TEtranscripts` requires pairs (-t/-c) and does not fit a per-sample workflow |
| TE GTF | `te_locus.exon.gtf` | Must contain only exon lines and must include `family_id`/`class_id` |

Note: locus-level and subfamily-level values **are not exactly equal** (multi-mapping is distributed more finely); the difference must be noted in the report.

## Telescope

| Parameter | Value | Description |
|---|---|---|
| `--attribute` | `gene_id` | **Must be specified explicitly**; the default is `locus`, and the `locus` attribute of the unified GTF equals `gene_id`, so results are identical either way |
| `--max_iter` | 200 | Consistent with the 2019 paper |
| `--theta_prior` | 200000 | Default value, also used by the old full rerun |
| `--ncpu` | **1** | Must be 1: >1 uses the parallel loader, which calls `fetch(region)` on the BAM and needs an index; collated BAMs cannot be indexed (`Chromosome blocks not continuous`). Parallelism comes from "concurrent samples", not multi-threading within one sample |
| `--reassign_mode` | default exclude | Consistent with the old full rerun, for easy comparison |

Prerequisite: the BAM must be `samtools collate`d (same read pair adjacent), otherwise assign fails.
In samtools 1.18, `-n` means "number of temporary files", not "sort by name".

## SQuIRE

| Parameter | Value | Description |
|---|---|---|
| `Map --read_length` | 36 | Consistent with historical runs; must be changed for other data |
| `Count --em` | auto | Consistent with historical runs |
| `clean_folder` | `squire_clean` | Native Build output, not replaced |

## Resources

| Parameter | Value | Basis |
|---|---|---|
| `max_concurrent_telescope` | 24 | Back-calculated from measured `/dev/shm` and memory usage; median single-sample output 6.7 GB |
| `telescope_mem_mb` | 120000 | 720k-locus pilot measured about 4.8 GB/sample, about 120 GB at 24 concurrent |
| `max_concurrent_squire` | 6 | 10 threads × 6 = 60, within the single-machine limit |
