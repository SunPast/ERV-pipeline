# -----------------------------------------------------------------------------
# TEtranscripts -- single-sample counting uses TEcount (from the same package, per sample, no control needed)
#
# Key point: TEtranscripts itself can do locus level, provided the GTF has "one locus, one gene".
#   This is exactly what the Telescope paper means by "modifications to the annotation were made to
#   allow locus-specific quantification"; no line of tool code needs to change.
# -----------------------------------------------------------------------------

rule tetranscripts_locus:
    input:
        bam = lambda wc: bam_of(wc),
        gene_gtf = P["gene_gtf"],
        te_gtf = f"{OUTDIR}/annotation/te_locus.exon.gtf",
    output:
        f"{OUTDIR}/tetranscripts/{{cohort}}/{{sample}}/{{sample}}.cntTable",
    threads: 1
    resources:
        mem_mb = config["resources"]["tetranscripts_mem_mb"],
        # TEcount writes its sorted temp BAM (~3.5 GB) into the process cwd,
        # so the job runs from the temp root instead of /data3.
        shm_gb = tmp_shm_gb,
        hdd_gb = tmp_hdd_gb,
        hdd2_gb = tmp_hdd2_gb,
    params:
        outdir = lambda wc: sdir("tetranscripts", wc.cohort, wc.sample),
        tmp = tmp_dir,
        bam_root = P["bam_root"],
        extra_bind = " ".join(f"--bind {d}:{d}:ro" for d in P.get("bind_paths", [])),
        ann_dir = f"{OUTDIR}/annotation",
        ref_dir = os.path.dirname(P["gene_gtf"]),
        apptainer = apptainer(),
        sif = P["tetranscripts_sif"],
    log:
        f"{OUTDIR}/tetranscripts/{{cohort}}/{{sample}}/run.log",
    shell:
        """
        tmp={params.tmp}.$$ && mkdir -p {params.outdir} "$tmp" && cd "$tmp" && \
        export TMPDIR="$tmp" && \
        {params.apptainer} exec \
            --bind {params.bam_root}:{params.bam_root}:ro \
            {params.extra_bind} \
            --bind {params.ann_dir}:{params.ann_dir}:ro \
            --bind {params.ref_dir}:{params.ref_dir}:ro \
            --bind {params.outdir}:{params.outdir} \
            --bind "$tmp":"$tmp" \
            {params.sif} TEcount --sortByPos --format BAM --mode multi \
            -b {input.bam} \
            --GTF {input.gene_gtf} --TE {input.te_gtf} \
            --project {wildcards.sample} --outdir {params.outdir} >> {log} 2>&1 && \
        cd / && rm -rf "$tmp"
        """
