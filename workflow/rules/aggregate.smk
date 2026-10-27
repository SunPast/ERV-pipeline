# -----------------------------------------------------------------------------
# Aggregation layer: four tools -> unified long table + subfamily evidence table
#
# Why not compute on the frontend:
#   The unified long table has tens of millions of rows, so the frontend cannot aggregate it on the fly; the evidence table
#   also uses cross-tool definitions and must be fixed in the backend, otherwise the numbers change on every refresh.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Standard tables consumed by the database loader
#
# The SQL side (ERV-data/rerun_out/db_load/build_import_csvs.sh) reads one CSV
# per tool/level, not the unified long table. This rule emits those files with
# unified names so only the paths in the loader need updating, not its logic.
# -----------------------------------------------------------------------------
rule export_standard:
    input:
        samples = SAMPLES_FILE,
        locus_map = f"{OUTDIR}/annotation/te_locus_map.tsv",
        ervmap_map = f"{OUTDIR}/annotation/ervmap_provirus_map.tsv",
        long = f"{OUTDIR}/aggregate/erv_expression_long.tsv.gz",
    output:
        directory(f"{OUTDIR}/standard"),
    params:
        outdir = OUTDIR,
        stddir = f"{OUTDIR}/standard",
        script = "workflow/scripts/export_standard.py",
    log:
        f"{OUTDIR}/aggregate/export_standard.log",
    shell:
        "mkdir -p {params.stddir} && python {params.script} {params.outdir} "
        "{input.samples} {input.locus_map} {input.ervmap_map} {params.stddir} 2> {log}"

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
        qc = f"{OUTDIR}/aggregate/qc_per_sample.tsv",
    params:
        outdir = OUTDIR,
        min_count = config["filter"]["min_count"],
        script = "workflow/scripts/aggregate_expression.py",
    log:
        f"{OUTDIR}/aggregate/aggregate.log",
    shell:
        "python {params.script} {params.outdir} {input.samples} "
        "{input.locus_map} {input.ervmap_map} {params.min_count} 2> {log}"
