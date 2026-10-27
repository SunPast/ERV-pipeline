# -----------------------------------------------------------------------------
# Telescope -- full-annotation rerun (550 subfamilies / 720,871 fragments)
#
# Why rerun:
#   The old result had only 60 internal-region models (14,968 merged loci), the only hard gap among the four tools.
#   With the unified GTF it shares the same set of coordinates with TEtranscripts, leaving only algorithmic differences.
#
# Input requirement: the BAM must be collated (same read pair adjacent), otherwise assign fails.
# -----------------------------------------------------------------------------

rule telescope_assign:
    input:
        bam = lambda wc: bam_of(wc),
        gtf = f"{OUTDIR}/annotation/te_locus.gtf",
    output:
        report = f"{OUTDIR}/telescope/{{cohort}}/{{sample}}/{{sample}}.report.tsv",
    threads: 8   # used only for collate; assign itself is single-threaded (seq mode)
    resources:
        mem_mb = config["resources"]["telescope_mem_mb"],
    params:
        outdir = lambda wc: sdir("telescope", wc.cohort, wc.sample),
        max_iter = config["tools"]["telescope"]["max_iter"],
        ncpu = config["tools"]["telescope"]["ncpu"],
        theta_prior = config["tools"]["telescope"]["theta_prior"],
        env = run_env(),
    log:
        f"{OUTDIR}/telescope/{{cohort}}/{{sample}}/assign.log",
    shell:
        """
        mkdir -p {params.outdir} && cd {params.outdir} && \
        {params.env} \
        samtools collate -@ {threads} -n 8 -l 6 -o collated.bam {input.bam} && \
        telescope assign --outdir . --exp_tag {wildcards.sample} \
            --ncpu {params.ncpu} --max_iter {params.max_iter} --theta_prior {params.theta_prior} \
            collated.bam {input.gtf} >> {log} 2>&1 && \
        mv -f {wildcards.sample}-telescope_report.tsv {output.report} && \
        rm -f collated.bam
        """
