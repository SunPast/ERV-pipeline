# -----------------------------------------------------------------------------
# SQuIRE
#
# Decision (2026-10-27, verified on real data; see docs/annotation.md):
#   * SQuIRE needs the FULL repeat annotation as a reference set, not ERV only.
#     With ERV-only input the counts of high-copy subfamilies collapse
#     (MSTA 622 -> 57) because the multi-mapper bookkeeping loses its anchor.
#   * The ERV part of the annotation is identical to our canonical loci
#     (720,871 rows, 0 coordinate differences, verified by comm).
#   * `squire Count` runs directly on the pipeline STAR BAM. The `squire Map`
#     step is therefore optional and disabled by default: it needs ~32 GB per
#     sample (STAR index) and is the single largest cost in the project.
#     Verified: Count on a plain coordinate-sorted STAR BAM exits 0 and writes
#     both TEcounts and subFcounts.
#
# Outputs (both levels come from ONE Count run):
#   <sample>_TEcounts.txt    locus level
#   <sample>_subFcounts.txt  subfamily level (1,180 rows; keep the 550 ERV ones)
#
# Caveats:
#   * The map folder must NOT contain any *.log file. `getlibsize()` picks the
#     first "*.log" and parses STAR log lines from it; a stray log (e.g. our own
#     stdout redirect) makes it crash with
#     "UnboundLocalError: local variable 'unique_libsize' referenced before assignment".
#   * SQuIRE TE_IDs use the format
#       chr|start|end|Subfamily:Family:Class|milliDiv|strand
#     which maps one-to-one to canonical locus_id by coordinates.
# -----------------------------------------------------------------------------

rule squire_count:
    input:
        bam = lambda wc: bam_of(wc),
        # The clean folder must exist first; without these inputs the job can
        # start before squire_clean_folder and fail with a missing *.bed.
        clean_bed = f"{OUTDIR}/annotation/squire_clean/hg38_all.bed",
        clean_copies = f"{OUTDIR}/annotation/squire_clean/hg38_all_copies.txt",
    output:
        te   = f"{OUTDIR}/squire/{{cohort}}/{{sample}}/{{sample}}_TEcounts.txt",
        subf = f"{OUTDIR}/squire/{{cohort}}/{{sample}}/{{sample}}_subFcounts.txt",
    threads: 6
    resources:
        mem_mb = config["resources"]["squire_mem_mb"],
    params:
        outdir = lambda wc: sdir("squire", wc.cohort, wc.sample),
        ref = P["squire_ref"],
        # Built by rule squire_clean_folder from the full-TE map, on local disk.
        clean = f"{OUTDIR}/annotation/squire_clean",
        apptainer = apptainer(),
        sif = P["squire_sif"],
        bam_root = P["bam_root"],
        extra_bind = " ".join(f"--bind {d}:{d}:ro" for d in P.get("extra_bind", [])),
        ann_dir = f"{OUTDIR}/annotation",
    shell:
        """
        mkdir -p {params.outdir} && \
        ln -sfn {input.bam} {params.outdir}/{wildcards.sample}.bam && \
        {params.apptainer} exec \
            --bind {params.bam_root}:{params.bam_root}:ro \
            {params.extra_bind} \
            --bind {params.ref}:{params.ref}:ro \
            --bind {params.ann_dir}:{params.ann_dir}:ro \
            --bind {params.outdir}:{params.outdir} \
            {params.sif} squire Count \
            --map_folder {params.outdir} \
            --clean_folder {params.clean} \
            --count_folder {params.outdir} \
            --read_length 36 --fetch_folder {params.ref} \
            --pthreads {threads} --name {wildcards.sample}
        """


# -----------------------------------------------------------------------------
# Optional: run `squire Map` when a FASTQ-based BAM is required.
# Enable with tools.squire.map: true. Costs ~47 min and ~32 GB per sample.
# -----------------------------------------------------------------------------
rule squire_map:
    input:
        fq1 = lambda wc: fq_of(wc)[0],
        fq2 = lambda wc: fq_of(wc)[1],
    output:
        bam = f"{OUTDIR}/squire_map/{{cohort}}/{{sample}}/{{sample}}.bam",
    threads: 10
    resources:
        mem_mb = 34000,
    params:
        outdir = lambda wc: f"{OUTDIR}/squire_map/{wc.cohort}/{wc.sample}",
        ref = P["squire_ref"],
        apptainer = apptainer(),
        sif = P["squire_sif"],
        fq_root = P["fastq_root"],
    shell:
        """
        mkdir -p {params.outdir} && \
        {params.apptainer} exec \
            --bind {params.fq_root}:{params.fq_root}:ro \
            --bind {params.ref}:{params.ref} \
            --bind {params.outdir}:{params.outdir} \
            {params.sif} squire Map \
            --read1 {input.fq1} --read2 {input.fq2} \
            --map_folder {params.outdir} --read_length 36 \
            --fetch_folder {params.ref} --pthreads {threads} \
            --name {wildcards.sample}
        """
