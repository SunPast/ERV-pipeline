# Deployment notes

Everything here is **host-specific**. A standard machine needs none of it — the
minimal config plus a `.sif` container is enough. These notes record what had to
be done on our node and why, so the same problems are not rediscovered.

## Container images: `.sif` vs sandbox directory

The pipeline expects `squire_sif` / `tetranscripts_sif`. On most hosts a `.sif`
image works directly. On a host **without squashfuse**, every `apptainer exec
<sif>` first unpacks the whole image (~2 GB) into `$TMPDIR`; with 160 concurrent
jobs that filled a 70 GB root filesystem within minutes and killed the run.

Workaround used here: build the image once as a sandbox directory and point the
config at the directory instead of the file.

```bash
apptainer build --sandbox /path/sandbox/squire  squire.sif
apptainer build --sandbox /path/sandbox/tetranscripts  tetranscripts.sif
```

## Scratch space

`tmp_roots` is optional. Empty (`[]`) means everything goes to `{outdir}/_tmp`,
which is correct but slow. On a multi-disk node, list the roots with shares; each
root is a separate Snakemake resource, so the job count per disk stays bounded.

Two measurements that shaped the defaults:

| Observation | Consequence |
|---|---|
| All scratch on one HDD, 160-way concurrency → 46 MB/s, 60/17,000 steps in 6 h | Never put every job's scratch on one spinning disk |
| 70 % of scratch moved to HDDs, 171 concurrent → throughput halved (22 vs 44 jobs/h) | RAM (`/dev/shm`) must stay the primary tier |

Scratch per job scales with input size: `3.5 × BAM`, floor 6 GB. A flat budget
under-sized the peak (16.6 GB observed for a 5 GB BAM) and tmpfs filled up,
which made intermediate files disappear mid-run.

## NFS input

BAMs usually live on shared storage. Measured: one stream 376 MB/s, eight
streams 1.4 GB/s, while the whole pipeline only pulled ~110 MB/s — the network
was never the bottleneck, and staging BAMs to local disk did **not** speed up
individual jobs. Do not build a staging layer on the assumption that it will.

What does show up: CPU sits ~77 % idle with `iowait ≈ 0` while jobs wait on NFS
replies, so the machine is limited by per-job blocking, not by core count. Adding
jobs beyond ~130 lowers total throughput.

## Resource caps

Passed at launch, not baked into the Snakefile, because they depend on the host:

```bash
snakemake -p --cores 160 \
  --resources shm_gb=480 hdd_gb=90 hdd2_gb=90 \
  --keep-going --rerun-incomplete all
```

`shm_gb` must stay below the actual size of `/dev/shm`; exceeding it lets jobs
reserve more than the tmpfs can hold.
