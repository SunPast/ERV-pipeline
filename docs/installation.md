# Installation

## Requirements

| Component | Version used here | Notes |
|---|---|---|
| Linux | RHEL 9 (x86_64) | needs `/dev/shm` larger than the largest per-sample temp file |
| Python | 3.14 (virtualenv in `.venv/`) | `snakemake >= 8`, `pandas` |
| Snakemake | 9.27.1 | `--keep-going`, resource counters |
| Apptainer | 1.x | only for SQuIRE and TEtranscripts images |
| samtools | 1.18 | collate for Telescope, BAM inspection |
| subread / featureCounts | 2.0.x | gene-level expression stage |
| Telescope | 1.0.3 | Python install under `paths.telescope_bin` |
| SQuIRE | 0.9.9.92 | inside `paths.squire_sif` |
| TEtranscripts | 2.2.3 | inside `paths.tetranscripts_sif` |

Python side:

```bash
python -m venv .venv
.venv/bin/pip install snakemake pandas
```

## Reference data

All paths live under `paths:` in `config/config.yaml`. Nothing is downloaded at run time.

| Config key | Content |
|---|---|
| `genome_fa` | GRCh38 primary assembly FASTA, with `.fai` |
| `gene_gtf` | GENCODE v34 annotation (gene level, for featureCounts) |
| `rmsk_table` | UCSC `hg38_rmsk.txt`; the single source of ERV coordinates and `milliDiv` |
| `te_gtf_official` | TEtranscripts hg38 RepeatMasker GTF (full TE, for SQuIRE) |
| `te_gtf_locus`, `te_exon_gtf` | ERV-only annotation, built by `rules/annotation.smk` |
| `ervmap_bed` | ERVmap's official 3,220 provirus BED |

The ERV annotation is not shipped: it is derived from `rmsk_table` on first run
(`annotation/te_locus.gtf`, 693,078 loci) and audited against the SQuIRE subset.

## Container images

`apptainer exec <file.sif>` unpacks the whole image into `$TMPDIR` unless
`squashfuse` is installed. On a host without it, 100+ concurrent jobs each extract
~2 GB and fill the root filesystem. Build a sandbox once and point the config at the
directory:

```bash
apptainer build --sandbox /path/to/sandbox/squire       /path/to/squire.sif
apptainer build --sandbox /path/to/sandbox/tetranscripts /path/to/tetranscripts.sif
```

```yaml
paths:
  squire_sif:        /path/to/sandbox/squire
  tetranscripts_sif: /path/to/sandbox/tetranscripts
```

## Configure

```bash
cp config/samples.example.tsv config/samples.tsv
$EDITOR config/config.yaml     # paths, outdir, tmp_roots, resources
$EDITOR config/samples.tsv     # cohort, sample, run, bam, fq1, fq2
```

`outdir` must be on node-local disk, not on the shared NFS mount. `tmp_roots`
lists the scratch directories and the share of jobs that use each; see
[usage.md](usage.md).

## Verify

```bash
.venv/bin/snakemake -n -p --cores 8 all
```

The dry run must report a job count on the order of `3 x samples` for the
quantification tools, plus one gene-count job per sample. A much smaller number
means the sample sheet was not picked up.
