# -----------------------------------------------------------------------------
# Gene-level expression
#
# Counts come from the same coordinate-sorted BAMs the ERV tools use, with the
# parameters of the TCASIA pipeline (workflows/01_alignment, rule featurecounts):
#
#   featureCounts -T {threads} -p -t exon -g gene_id -a <gencode v34> -o <out> <bam>
#
# -s is left at 0 (unstranded) on purpose: the existing TCASIA matrices were
# built that way, so new samples stay comparable with them.
#
# No expression filter is applied here. Filtering is an analysis-time choice
# (DESeq2 independent filtering, edgeR filterByExpr) and baking one rule into
# the platform would make other analyses unreproducible.
# -----------------------------------------------------------------------------

rule gene_counts:
    input:
        bam = bam_of,
        gtf = P["gene_gtf"],
    output:
        f"{OUTDIR}/gene_counts/{{cohort}}/{{sample}}/{{sample}}_featureCounts.txt",
    threads: config["tools"]["gene_counts"]["threads"]
    params:
        bin = P["featurecounts_bin"],
        tmp = lambda w: f"{tmp_dir(w)}_genecounts",
    resources:
        shm_gb  = tmp_shm_gb,
        hdd_gb  = tmp_hdd_gb,
        hdd2_gb = tmp_hdd2_gb,
    log:
        f"{OUTDIR}/gene_counts/logs/{{cohort}}/{{sample}}.log",
    shell:
        """
        mkdir -p {params.tmp} $(dirname {log:q})
        export TMPDIR={params.tmp}
        {params.bin}/featureCounts -T {threads} -p -t exon -g gene_id \
            -a {input.gtf:q} -o {output:q} {input.bam:q} >> {log} 2>&1
        rm -rf {params.tmp}
        """


rule merge_gene_counts:
    input:
        samples = SAMPLES_FILE,
        gtf = P["gene_gtf"],
        counts = gene_count_outputs(),
    output:
        std = f"{OUTDIR}/gene_expression/gene_expression_standard.csv.gz",
        matrices = directory(f"{OUTDIR}/gene_expression/matrices"),
        factors = f"{OUTDIR}/gene_expression/sample_gene_factor.tsv",
    params:
        outdir = OUTDIR,
        script = "workflow/scripts/merge_gene_counts.py",
        factor_script = "workflow/scripts/compute_gene_factors.py",
    log:
        f"{OUTDIR}/gene_expression/merge.log",
    shell:
        """
        python {params.script} {params.outdir} {input.samples} {input.gtf} 2>> {log}
        python {params.factor_script} {params.outdir} {input.samples} 2>> {log}
        """
