# Usage

## Sample sheet

Tab separated, one row per sample. `cohort` is the label used everywhere
downstream; it does not have to match the upstream directory name.

```text
cohort	sample	run	bam	fq1	fq2
PHS002176	SRR13780250	SRR13780250	/data/data2/IO_RNA/AS/PHS002176/aligned/SRR13780250/SRR13780250_Aligned.sortedByCoord.out.bam		
```

`bam` must exist and be non-empty; `make_samples.py` skips zero-byte files and
falls back to a search when a path is missing. `fq1`/`fq2` are only needed when
`squire.map` is enabled.

## Running

```bash
.venv/bin/snakemake -n -p --cores 8 all          # dry run
.venv/bin/snakemake -p --cores 160 --keep-going all
```

`--keep-going` keeps the run alive when individual samples fail (for example a
BAM with no uniquely mapped reads); failures are counted at the end and can be
retried by re-running the same command, `--rerun-incomplete` picks up
half-written outputs.

Targets:

| Target | Content |
|---|---|
| `all` | every per-sample output + aggregate + standard tables + gene expression |
| `annotation/te_locus.gtf` | ERV annotation only, no quantification |
| `standard` | per-tool CSVs for the database loader |
| `gene_expression/gene_expression_standard.csv.gz` | gene-level table only |

## Scratch space

Each of the four tools writes GB-scale temporary files. `tmp_roots` assigns a
sample to one root by a hash of its name, so a resumed run keeps the same
placement, and `resources.tmp_gb_per_job` charges the job against that root's
counter.

```yaml
tmp_roots:
  - {path: /dev/shm/ervtmp,      share: 0.80, resource: shm_gb}
  - {path: /data4/.../ervtmp,    share: 0.10, resource: hdd_gb}
  - {path: /data3/.../ervtmp,    share: 0.10, resource: hdd2_gb}
```

Launch with matching counters, e.g. `--resources shm_gb=450 hdd_gb=300 hdd2_gb=300`.
Putting all jobs on one HDD collapses it to tens of MB/s at high concurrency.

## Outputs

| Path | Content |
|---|---|
| `telescope/<cohort>/<sample>/<sample>.report.tsv` | per-locus counts |
| `squire/<cohort>/<sample>/<sample>_{TE,subF}counts.txt` | locus and subfamily counts |
| `tetranscripts/<cohort>/<sample>/<sample>.cntTable` | locus counts |
| `ervmap/<cohort>/<sample>/<sample>.locus.tsv` | imported ERVmap result |
| `gene_counts/<cohort>/<sample>/<sample>_featureCounts.txt` | gene-level counts |
| `aggregate/erv_expression_long.tsv.gz` | unified ERV long table |
| `aggregate/erv_subfamily_evidence.tsv.gz` | cross-tool subfamily evidence |
| `aggregate/qc_per_sample.tsv` | counts and positive loci per sample per tool |
| `standard/*.csv` | inputs of the database loader |
| `gene_expression/gene_expression_standard.csv.gz` | gene-level table for the loader |

Field definitions: [output_format.md](output_format.md).

## Checking a run

```bash
grep -c '^Finished jobid' logs/full_*.log     # completed steps
grep '^Error in rule' logs/full_*.log | sort | uniq -c
column -t aggregate/qc_per_sample.tsv | awk '$NF!="OK"'
```

`qc_per_sample.tsv` flags `ZERO_COUNTS` and `LOW_COUNTS`; a run that finishes with
errors should be re-invoked with the same command, completed samples are not
recomputed.
