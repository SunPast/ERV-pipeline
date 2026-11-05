# ERV-pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-9.27.1-1f425f.svg)](https://snakemake.readthedocs.io)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Quantify human endogenous retroviruses (ERVs) and transposable elements from
RNA-seq alignment files with four independent tools, then put their output on a
common annotation so the results can be compared.

| Tool | Level produced | Annotation it reads |
|---|---|---|
| [TEtranscripts](https://github.com/mhammell-laboratory/TEtranscripts) | locus + subfamily | ERV exon GTF (693,078 loci, 550 subfamilies) |
| [Telescope](https://github.com/mlbendall/telescope) | locus + subfamily | same 693,078 loci, gene lines only |
| [SQuIRE](https://github.com/wyang1989/SQuIRE) | locus + subfamily | full TE annotation (4,520,928 loci) |
| [ERVmap](https://github.com/mtokuyama/ERVmap) | provirus + subfamily | the official 3,220-provirus BED |

The four tools are deliberately given the annotation form each one requires. ERVmap's
BED is the method itself, so it is imported rather than recomputed; the other three
are quantified here from a coordinate-sorted BAM.

Gene-level expression is produced from the same BAMs with `featureCounts`, so ERV
and gene counts always describe the same library.

## Requirements

- Linux with a container runtime (`apptainer` or `singularity`) for the SQuIRE and
  TEtranscripts images
- Snakemake ≥ 9
- Python ≥ 3.10 with `pandas`
- Reference files: genome FASTA (+ `.fai`), a GENCODE GTF, the UCSC RepeatMasker
  table, and the public TEtranscripts hg38 GTF
- Input: one coordinate-sorted BAM per sample

Exact versions and the reference layout are in [docs/installation.md](docs/installation.md).

## Quick start

```bash
cp config/config.minimal.yaml config/config.yaml   # set refdir and outdir
cp config/samples.example.tsv config/samples.tsv   # cohort, sample, bam per line
snakemake -n -p                                    # dry-run: check the plan
snakemake -p --cores 32 all                        # run
```

Reference paths are derived from `refdir`; list a path under `paths:` only when a
file does not follow the standard layout. Everything host-specific (container
sandboxes, scratch layout, resource caps) is in
[docs/deployment.md](docs/deployment.md) and is not needed on a standard machine.

Run a single layer instead of everything:

```bash
snakemake -p --cores 8 annotation/te_locus.gtf annotation/annotation_check.tsv
```

## Input modes

The pipeline starts from either FASTQ or BAM; no external alignment step is required.

| samples.tsv | What runs |
|---|---|
| `fq1` + `fq2` filled | `fastp` → `STAR` → `samtools sort/index` → four quantification tools |
| `bam` filled | quantification tools only (legacy cohorts) |

Alignment is paired-end only. Every supplied BAM is validated before use:
`*.tmp.*` STAR shards are rejected, and a missing BGZF EOF block (truncated file)
aborts the run instead of silently producing partial counts.

## Configuration

`config/config.yaml` holds paths, per-tool switches and resource limits. Each tool
can be turned off independently (`tools.<name>.run`), and
`paths.external_gene_dir` optionally brings in pre-computed gene matrices instead
of quantifying those samples again. See [docs/parameters.md](docs/parameters.md).

## Output

- `aggregate/erv_expression_long.tsv.gz` — one long table, all tools and both
  levels: `tool, cohort, sample, level, feature, region_type, count, milli_div`
  (`milli_div` is filled at locus level only)
- `aggregate/erv_subfamily_evidence.tsv.gz` — per-subfamily support across tools
- `standard/` — per-tool CSV files ready for a database loader
- `gene_expression/` — long table for the loader plus one matrix per cohort
- `annotation/annotation_check.tsv` — differences between the annotations used

Field definitions are in [docs/output_format.md](docs/output_format.md).

## Documentation

| File | Contents |
|---|---|
| [docs/installation.md](docs/installation.md) | dependencies, reference data, container images |
| [docs/usage.md](docs/usage.md) | sample sheet, running and resuming, outputs |
| [docs/annotation.md](docs/annotation.md) | why each tool gets a different annotation, with measurements |
| [docs/literature.md](docs/literature.md) | how other pipelines handle annotation |
| [docs/gene_expression.md](docs/gene_expression.md) | gene counts, length definition, optional import mode |
| [docs/parameters.md](docs/parameters.md) | parameter reference |
| [docs/troubleshooting.md](docs/troubleshooting.md) | pitfalls established by measurement |
| [docs/deployment.md](docs/deployment.md) | host-specific notes |
| [CHECKLIST.md](CHECKLIST.md) | pre-run and post-run checks |

## Scope

Subfamily level is comparable across all four tools. Locus level is mapped onto a
common coordinate set but never re-allocated between loci; differences between
tools at that level are expected and should be read as such.

## Related

- [ERV-motif](https://github.com/SunPast/ERV-motif) — static annotation of
  transcription-factor binding motifs in each ERV subfamily's consensus
  sequence. It has no sample dimension and joins the tables produced here on
  `subfamily`.

## Citation

Cite the four tools and the reference releases listed in [CITATION.md](CITATION.md).

## License

MIT — see [LICENSE](LICENSE).

### Alignment parameters

The FASTQ stage mirrors TCASIA_pipeline (`workflows/01_alignment`) so that BAMs and
gene counts are directly comparable:

| Step | Parameters |
|---|---|
| fastp | `--length_required 36 --qualified_quality_phred 20 --unqualified_percent_limit 40 --n_base_limit 5`, `--detect_adapter_for_pe` for paired-end |
| STAR | `--outSAMunmapped Within --outFilterMismatchNmax 15 --outSAMattributes All --outFilterIntronMotifs RemoveNoncanonical --outSAMstrandField intronMotif --twopassMode Basic --quantMode GeneCounts --outSAMtype BAM Unsorted` |
| samtools | `sort` then `index` |
| featureCounts | `-p -t exon -g gene_id` (unstranded) |

TCASIA_pipeline requires paired-end input. Single-end is supported here as an
extension: `--detect_adapter_for_pe` is dropped and STAR receives one read file.

