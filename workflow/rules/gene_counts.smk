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
        # -O: allow a read to be counted for every gene it overlaps.
        # Without it, featureCounts drops shared reads as Unassigned_Ambiguity;
        # EEF1G is fully covered by the lncRNA AP002990.1 on the same strand and
        # ended up with 1-6 counts in 171/4289 samples (69,892 reads in the BAM).
        {params.bin}/featureCounts -T {threads} -p -O -t exon -g gene_id \
            -a {input.gtf:q} -o {output:q} {input.bam:q} >> {log} 2>&1
        rm -rf {params.tmp}
        """


rule merge_gene_counts:
    input:
        samples = SAMPLES_FILE,
        gtf = P["gene_gtf"],
        counts = gene_count_outputs(),
        # empty unless paths.external_gene_dir is set
        external = external_gene_outputs(),
    output:
        std = f"{OUTDIR}/gene_expression/gene_expression_standard.csv.gz",
        matrices = directory(f"{OUTDIR}/gene_expression/matrices"),
        factors = f"{OUTDIR}/gene_expression/sample_gene_factor.tsv",
    params:
        outdir = OUTDIR,
        external = (external_gene_outputs() or [""])[0],
        script = "workflow/scripts/merge_gene_counts.py",
        factor_script = "workflow/scripts/compute_gene_factors.py",
    log:
        f"{OUTDIR}/gene_expression/merge.log",
    shell:
        """
        python {params.script} {params.outdir} {input.samples} {input.gtf} {params.external} 2>> {log}
        python {params.factor_script} {params.outdir} {input.samples} 2>> {log}
        """
