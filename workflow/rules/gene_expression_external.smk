# -----------------------------------------------------------------------------
# External gene-count matrices
#
# Only included when paths.external_gene_dir is set. Some cohorts already have
# per-cohort count matrices computed outside this repository; recomputing them
# costs hours and adds nothing. This rule converts them into the pipeline's long
# format and keeps a `source` column, so a reader can always tell which numbers
# were computed here and which were imported.
# -----------------------------------------------------------------------------
rule import_external_gene_counts:
    """Convert pre-computed cohort matrices into the pipeline's long format."""
    output:
        tsv = f"{OUTDIR}/gene_expression/external_counts.tsv.gz",
    params:
        dir     = P["external_gene_dir"],
        pattern = config.get("external_gene_pattern", "{cohort}_merged_counts.csv"),
        samples = SAMPLES_FILE,
        source  = config.get("external_gene_source", "external"),
        script  = "workflow/scripts/import_external_gene_counts.py",
    log:
        f"{OUTDIR}/gene_expression/logs/import_external.log",
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log:q})
        rm -f {output.tsv}
        python3 {params.script} \
            --dir {params.dir:q} \
            --pattern '{params.pattern}' \
            --samples {params.samples:q} \
            --out {output.tsv} \
            --source-name {params.source} 2> {log}
        """
