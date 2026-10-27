# -----------------------------------------------------------------------------
# Aggregation layer: four tools -> unified long table + subfamily evidence table
#
# Why not compute on the frontend:
#   The unified long table has tens of millions of rows, so the frontend cannot aggregate it on the fly; the evidence table
#   also uses cross-tool definitions and must be fixed in the backend, otherwise the numbers change on every refresh.
# -----------------------------------------------------------------------------

rule aggregate_expression:
    input:
        samples = SAMPLES_FILE,
        locus_map = f"{OUTDIR}/annotation/te_locus_map.tsv",
        ervmap_map = f"{OUTDIR}/annotation/ervmap_provirus_map.tsv",
        # Per-sample outputs of the rerun tools. Listing them here is what makes
        # `snakemake all` actually execute the quantification rules; aggregation
        # itself only reads whatever is on disk.
        tool_outputs = all_tool_outputs(),
    output:
        long = f"{OUTDIR}/aggregate/erv_expression_long.tsv.gz",
        evid = f"{OUTDIR}/aggregate/erv_subfamily_evidence.tsv.gz",
    params:
        outdir = OUTDIR,
        min_count = config["filter"]["min_count"],
        script = "workflow/scripts/aggregate_expression.py",
    log:
        f"{OUTDIR}/aggregate/aggregate.log",
    shell:
        "python {params.script} {params.outdir} {input.samples} "
        "{input.locus_map} {input.ervmap_map} {params.min_count} 2> {log}"
