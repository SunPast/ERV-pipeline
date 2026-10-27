# ERV-pipeline

A unified processing pipeline for four ERV/TE quantification tools (ERVmap, Telescope, SQuIRE, TEtranscripts).
Structure follows [TCASIA_pipeline](https://github.com/OncoHarmony-Network/TCASIA_pipeline): modular rules + environment declarations + documentation + validation report.

## Three ground rules

1. **One reference universe**: hg38 + GENCODE v34 + one UCSC rmsk table. Every ERV locus, its subfamily and its milliDiv come from the same release.
2. **Inputs are not uniform**: each tool gets the form it accepts; we do not force the same file onto all of them.
3. **Comparison is layered**: subfamily level is comparable; locus level is only mapped, never numerically redistributed.

## How the four tools are handled

| Tool | Annotation source | Rerun | Granularity | Notes |
|---|---|---|---|---|
| **TEtranscripts** | `te_locus.exon.gtf` (693,078 ERV loci) | **Yes** | locus + subfamily | `TEcount` per sample; subfamily aggregated from locus |
| **Telescope** | `te_locus.gtf` (same 693,078 loci, gene lines only) | **Yes** | locus + subfamily | Old result covered only 60 internal-region models |
| **SQuIRE** | `annotation/squire_clean` (full TE, 4,520,928 loci) | **Yes** | locus + subfamily | Built from the same rmsk table, not from `squire Build`; one `Count` run writes both levels |
| **ERVmap** | Official 3,220 provirus BED | No | locus (provirus) + subfamily (mapped) | The BED is the method itself; changing the annotation invalidates it |

Why SQuIRE needs the **full** TE annotation (and not just the ERV subset) is measured in
[docs/annotation.md](docs/annotation.md): with ERV-only input, counts of high-copy subfamilies collapse
(e.g. MSTA 622 -> 57) because its multi-mapper bookkeeping loses its reference set.

`SQuIRE Map` is **not** run: `squire Count` reads the existing coordinate-sorted STAR BAM directly.
That removes the 32 GB/sample STAR index from the picture.

## Directory layout

```
config/
  config.yaml          all paths/parameters/tool strategies
  samples.example.tsv  sample sheet template (cohort, sample, run, bam, fq1, fq2)
workflow/
  Snakefile
  rules/
    common.smk         shared variables and functions
    annotation.smk     annotation build from the rmsk table + SQuIRE clean folder + audit
    tetranscripts.smk  locus-level quantification
    telescope.smk      collate + assign
    squire.smk         Count on the STAR BAM (both levels)
    ervmap.smk         result import + provirus composition mapping
    aggregate.smk      aggregate into long table + evidence table
  scripts/             pure Python, runnable standalone
docs/
  annotation.md        annotation sources, decisions and measurements (most important)
  literature.md        what other pipelines do about annotation
  output_format.md     output table structure
  parameters.md        parameter reference
  troubleshooting.md   empirically established pitfalls
  work_order.md        full-rerun work order (node/storage policy, phases, acceptance)
```

## Quick start

```bash
cd ~/ERV-pipeline
cp config/samples.example.tsv config/samples.tsv    # fill in with real samples
snakemake -n -p                                     # dry-run
snakemake -p --cores 160 all
```

`outdir` in `config/config.yaml` points at node-local disk; products are not written to the
master NFS mount. See `docs/work_order.md` for the storage/migration policy.

Run only the annotation layer:

```bash
snakemake -p --cores 8 \
  annotation/te_locus.gtf \
  annotation/annotation_check.tsv \
  annotation/ervmap_provirus_map.tsv
```

## Output

- `aggregate/erv_expression_long.tsv.gz`: unified long table
  (tool, cohort, sample, level, feature, region_type, count, milli_div)
  `milli_div` is filled at **locus level only**; subfamily and ERVmap rows leave it empty.
- `aggregate/erv_subfamily_evidence.tsv.gz`: subfamily evidence table (n_tools_support, ervmap_support)
- `standard/`: six per-tool CSV files for the database loader (see below)
- `annotation/annotation_check.tsv`: difference metrics across annotations

Field definitions are in [docs/output_format.md](docs/output_format.md).

## Relationship to the online database

The SQL side does **not** read the long table. `ERV-data/rerun_out/db_load/build_import_csvs.sh`
reads one CSV per tool/level, then maps sample names to pseudonyms, adds `norm_count`
from its own size factors and loads via `load.sh`. This pipeline therefore also emits
those files (`standard/`, unified names):

| Standard file | Columns | Loads into |
|---|---|---|
| `ervmap_standard.csv` | cohort, sample, chr, start, end, ervmap_id, strand, count | `tool_ervmap` |
| `tetranscripts_subf_standard.csv` | cohort, sample, subfamily, family, count | `tool_tetranscripts` |
| `tetranscripts_locus_standard.csv` | cohort, sample, locus, subfamily, family, chr, start, end, strand, count | `tool_tetranscripts_locus` (new) |
| `telescope_locus_standard.csv` | same as above | `tool_telescope` |
| `squire_te_standard.csv` | cohort, sample, subfamily, family, count, TE_chr, TE_start, TE_stop, milliDiv, TE_strand | `tool_squire_te` |
| `squire_subf_standard.csv` | cohort, sample, subfamily, family, count | `tool_squire_subf` |

Two consequences for the SQL side:

1. `build_import_csvs.sh` must point at the new file names (`standard/` paths);
   the conversion logic itself is unchanged.
2. The locus-level TEtranscripts output needs a new table `tool_tetranscripts_locus`
   with the same columns as `tool_telescope`.

The loader keeps only families `ERV1 / ERVK / ERVL / ERVL-MaLR`, so the 3,387 loci
annotated as `ERVL?` or `ERV1?` (0.49%) are dropped at load time by design.

## Validation status (measured on this machine)

| Item | Result |
|---|---|
| ERV annotation (chr1-22/X/Y, one rmsk fragment per locus) | **693,078 loci**, 550 subfamilies |
| Full-TE annotation (SQuIRE reference) | **4,520,928 loci**, 1,180 subfamilies |
| Our ERV loci vs SQuIRE native ERV loci | 693,078 vs 693,078, **0 differences** (row-wise) |
| Annotation build time | 12 s |
| ERVmap composition mapping | 3,220 loci: 2,821 hit a provirus, 399 fell back to rmsk |
| Telescope, one sample, 15.7M fragments | **16 min 42 s**, peak **2.45 GB**, assignment rate 96.7% |
| Telescope output | 5,881 loci with records (5,194 with counts) |
| SQuIRE Count on a STAR BAM, 1M read pairs | **1 min 36 s**, peak **3.2 GB**, both output levels written |
| TEcount, one sample | **1 h 24 min**, peak **2.61 GB** |

## Key constraints (empirically established)

1. **Telescope must use `--ncpu 1`**: the parallel loader needs a BAM index, and collated BAMs cannot be indexed.
2. **Telescope's GTF must contain only gene lines**: it does not filter by feature type, and duplicated
   gene/transcript/exon entries change the result.
3. **TEcount requires `family_id`/`class_id`** in the TE GTF, otherwise it reports "TE GTF format error".
4. **SQuIRE Count accepts an external annotation** through `--clean_folder` (one `*.bed` + one `*_copies.txt`),
   but that folder must contain **no `*.log`** file: `getlibsize()` parses the first `*.log` it finds and
   crashes if the expected STAR lines are missing.
5. **Only primary chromosomes are kept** (chr1-22, X, Y). alt/fix/random contigs are named differently from
   the contigs in the BAM headers, so those annotation rows can never receive a read; chrM carries no ERV.
