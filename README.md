# ERV-pipeline

A unified processing pipeline for four ERV/TE quantification tools (ERVmap, Telescope, SQuIRE, TEtranscripts).
Structure follows [TCASIA_pipeline](https://github.com/OncoHarmony-Network/TCASIA_pipeline): modular rules + environment declarations + documentation + validation report.

## Three ground rules

1. **One reference universe**: hg38 + GENCODE v34 + a single rmsk (the official TEtranscripts hg38 RepeatMasker GTF).
2. **Inputs are not uniform**: each tool gets the form it accepts; we do not force the same file onto all of them.
3. **Comparison is layered**: subfamily level is comparable; locus level is only mapped, never numerically redistributed.

## How the four tools are handled

| Tool | Annotation source | Rerun | Granularity | Notes |
|---|---|---|---|---|
| **TEtranscripts** | Unified reference (rmsk fragment GTF) | **Yes** | locus + subfamily | Swap in a locus GTF to get locus level; subfamily is aggregated from locus |
| **Telescope** | Unified reference (same GTF) | **Yes** | locus + subfamily | Old result covered only 60 internal-region models, must be completed |
| **SQuIRE** | Native (`squire Build` auto-downloads UCSC rmsk) | No | locus + subfamily | Does not accept external annotation; median overlap with the reference is 96.8% |
| **ERVmap** | Official 3,220 provirus BED | No | locus (provirus) + subfamily (mapped) | The BED is the method itself; changing the annotation invalidates it |

For per-tool annotation differences, measured numbers, and "why they cannot be unified", see [docs/annotation.md](docs/annotation.md).

## Directory layout

```
config/
  config.yaml          all paths/parameters/tool strategies
  samples.example.tsv  sample sheet template (cohort, sample, run, bam, fq1, fq2)
workflow/
  Snakefile
  rules/
    common.smk         shared variables and functions
    annotation.smk     unified annotation build + validation of the three annotations
    tetranscripts.smk  locus-level quantification
    telescope.smk      collate + assign
    squire.smk         Map + Count (native annotation)
    ervmap.smk         result import + provirus composition mapping
    aggregate.smk      aggregate into long table + evidence table
  scripts/             pure Python, runnable standalone
docs/
  annotation.md        annotation sources and differences (most important)
  output_format.md     output table structure
  parameters.md        parameter reference
```

## Quick start

```bash
cd ~/ERV-pipeline
cp config/samples.example.tsv config/samples.tsv    # fill in with real samples
snakemake -n -p                                     # dry-run
snakemake -p --cores 32 all
```

Run only the annotation layer and aggregation (no quantification):

```bash
snakemake -p --cores 8 \
  results/annotation/te_locus.gtf \
  results/annotation/annotation_check.tsv \
  results/annotation/ervmap_provirus_map.tsv
```

## Output

- `results/aggregate/erv_expression_long.tsv.gz`: unified long table across the four tools (tool, cohort, sample, level, feature, region_type, count)
- `results/aggregate/erv_subfamily_evidence.tsv.gz`: subfamily evidence table (n_tools_support, ervmap_support)
- `results/annotation/annotation_check.tsv`: difference metrics across the three annotations

Field definitions are in [docs/output_format.md](docs/output_format.md).

## Relationship to the online database

This pipeline only produces files; it does not write to the database. The loading scripts remain in `ERV-data/rerun_out/db_load/sql/`;
the two are connected via `erv_expression_long.tsv.gz`, avoiding a situation where "the pipeline is one thing and the database loading another".

## Validation status (measured on this machine)

| Item | Result |
|---|---|
| Annotation layer build (720,871 loci) | 12.6 s, 550 subfamilies |
| Validation of the three annotations | reference 720,871 / SQuIRE 702,850; per-subfamily overlap median **96.8%**, 399 subfamilies ≥95% |
| ERVmap composition mapping | 3,220 loci: 2,821 hit a provirus, 399 fell back to rmsk; 2,646 have a 5'LTR, 2,066 have a 3'LTR |
| Telescope single sample (720k loci, single thread) | **15 min 52 s**, peak memory **2.57 GB**, 15.7M fragments |
| Telescope output | 5,881 loci with records (5,194 with counts); read assignment rate is on the same order as the old 60-model run |
| Snakemake full-pipeline dry-run | Passed; annotation + aggregation actual run 19 s |
| TEcount locus level (same sample) | **87 minutes**, peak **2.59 GB**; 720,871 output rows, 3,381 loci with counts, total TE counts 20,519, covering 390 subfamilies |

## Key constraints (empirically established)

1. **Telescope must use `--ncpu 1`**: the parallel loader needs a BAM index, and collated BAMs cannot be indexed.
2. **Telescope's GTF must contain only gene lines**: it does not filter by feature type, and duplicated gene/transcript/exon entries change the result (measured: gene-only 5,881 loci vs gene+transcript 4,045).
3. **TEcount requires `family_id`/`class_id`**: otherwise it reports "TE GTF format error".
4. **SQuIRE has no `MakeBed`**: 0.9.9.92 can only use `Build` to auto-download rmsk and does not accept external annotation.
