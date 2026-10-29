# Pre-run checklist

Everything below is checked on the host before a full run is started. Each item
corresponds to a failure that actually occurred during the 2026-10 rerun.

## Inputs

- [ ] `config/samples.tsv` lists every sample once; `wc -l` equals `samples + 1`
- [ ] every `bam` path exists and is **non-empty** (`find ... -size 0` returns nothing)
- [ ] BAM headers are coordinate sorted (`samtools view -H … | grep SO:coordinate`)
- [ ] contig names in the BAMs match the annotation (`chr1`, not `1`)
- [ ] all BAM roots are listed in `paths.bind_paths`; containers mount each root read-only

## Reference

- [ ] `rmsk_table`, `te_gtf_official`, `gene_gtf`, `genome_fa` exist and their
      SHA-256 values are recorded (`annotation_check.tsv` is regenerated)
- [ ] ERV annotation audits to 693,078 loci and 550 subfamilies
- [ ] the SQuIRE clean folder has no `*.log` file

## Runtime

- [ ] `outdir` is on node-local disk, with free space >= 3x the expected output size
- [ ] each `tmp_roots` entry exists and is writable
- [ ] `/dev/shm` is at least 400 GB, or `share` is reduced accordingly
- [ ] container images are pre-built sandboxes, not `.sif` files
- [ ] `--resources` counters passed at launch match `tmp_roots`
- [ ] no stale `.snakemake/locks` directory from an interrupted run

## Before the first job starts

- [ ] dry run reports the expected job count (`~3 x samples` + one gene-count job each)
- [ ] a single sample has been run end to end and its `qc_per_sample.tsv` row is `OK`
- [ ] `nproc` and `free -g` are recorded, and `--cores` leaves headroom for other users

## While running

- [ ] `grep '^Error in rule' log | sort | uniq -c` shows only known, sample-specific failures
- [ ] root filesystem stays below 80%; no `No space left on device`
- [ ] no tool writes to `/tmp` (GNU `sort` uses `$TMPDIR`)
- [ ] D-state process count stays low; rising `iowait` means scratch I/O contention

## After the run

- [ ] every sample has a row in `qc_per_sample.tsv`; `ZERO_COUNTS` rows are investigated
- [ ] standard tables and gene expression table are written
- [ ] products are rsync'ed to the archive location, then the database loader is run
