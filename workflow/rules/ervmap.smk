# -----------------------------------------------------------------------------
# ERVmap -- no rerun, import + mapping only
#
# Reason: its 3,220 provirus sequences are the method itself. After changing the annotation (e.g. to rmsk fragments),
#   the premise of "strict best-match + mismatch filtering" no longer holds, and the result is no longer ERVmap.
#
# Approach:
#   1) Put the existing ERVmap results as <cohort>/<sample>ERVresults.txt into results/ervmap_raw/
#      (scripts/link_ervmap_raw.sh can batch-symlink them from ERV-data/rerun)
#   2) This pipeline normalizes them into the unified long table and additionally produces a "provirus composition" mapping:
#      ERVmap locus -> provirus locus -> int_model / ltr5_model / ltr3_model / category
# -----------------------------------------------------------------------------

rule ervmap_import:
    input:
        raw = f"{OUTDIR}/ervmap_raw/{{cohort}}/{{sample}}ERVresults.txt",
        catalog = P["ervmap_catalog"],
    output:
        tsv = f"{OUTDIR}/ervmap/{{cohort}}/{{sample}}/{{sample}}.locus.tsv",
    params:
        script = "workflow/scripts/normalize_ervmap.py",
    log:
        f"{OUTDIR}/ervmap/{{cohort}}/{{sample}}/import.log",
    shell:
        "mkdir -p $(dirname {output.tsv}) && python {params.script} {input.raw} {input.catalog} {output.tsv} 2> {log}"

rule ervmap_provirus_map:
    """One-off: assignment table of the 3,220 provirus loci -> rmsk fragments / internal regions / flanking LTRs."""
    input:
        catalog = P["ervmap_catalog"],
        provirus = P["provirus_gtf"],
        lodict = f"{OUTDIR}/annotation/te_locus_map.tsv",
    output:
        f"{OUTDIR}/annotation/ervmap_provirus_map.tsv",
    params:
        script = "workflow/scripts/map_ervmap_provirus.py",
    log:
        f"{OUTDIR}/annotation/ervmap_provirus_map.log",
    shell:
        "mkdir -p $(dirname {output}) && python {params.script} {input.catalog} {input.provirus} {input.lodict} {output} 2> {log}"
