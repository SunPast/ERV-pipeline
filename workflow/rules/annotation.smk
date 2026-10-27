# -----------------------------------------------------------------------------
# Annotation layer
#
# Single source of truth: the official TEtranscripts hg38 RepeatMasker GTF.
# It is a FULL-TE file (4,693,511 loci / 1,180 subfamilies), not ERV only.
#
# Two derived products:
#   1) ERV subset (550 subfamilies, 720,871 loci) -> shared by Telescope and
#      TEtranscripts. One rmsk fragment = one locus, no merging.
#   2) Full-TE map (1,180 subfamilies) -> used to build the SQuIRE clean folder.
#      SQuIRE needs the complete repeat set as its multi-mapper reference;
#      with ERV-only input, counts of high-copy subfamilies collapse.
#
# Why fragments are NOT merged:
#   - merging rules are hand-made (distance threshold, model coverage) and
#     cannot be reproduced;
#   - fragments can always be aggregated upward, merged loci cannot be split
#     back into fragments;
#   - fragment-level coordinates are identical to TEtranscripts, so tool
#     differences can be attributed to the algorithm, not to the annotation.
# -----------------------------------------------------------------------------

rule build_locus_gtf:
    input:
        rmsk = P["rmsk_table"],
        dict = P["subfamily_dict"],
    output:
        gtf  = f"{OUTDIR}/annotation/te_locus.gtf",
        map  = f"{OUTDIR}/annotation/te_locus_map.tsv",
        exon = f"{OUTDIR}/annotation/te_locus.exon.gtf",
    log:
        f"{OUTDIR}/annotation/build_locus_gtf.log",
    params:
        script = "workflow/scripts/build_locus_from_rmsk.py",
    shell:
        "mkdir -p $(dirname {output.gtf}) && "
        "python {params.script} {input.rmsk} {input.dict} "
        "{output.gtf} {output.map} {output.exon} 2> {log}"


# Full-TE map only (no GTF needed): used by SQuIRE.
# Passing "-" as the dictionary disables the ERV filter.
rule build_full_locus_map:
    input:
        rmsk = P["rmsk_table"],
    output:
        map = f"{OUTDIR}/annotation/te_locus_map_full.tsv",
    log:
        f"{OUTDIR}/annotation/build_full_locus_map.log",
    params:
        script = "workflow/scripts/build_locus_from_rmsk.py",
    shell:
        "mkdir -p $(dirname {output.map}) && "
        "python {params.script} {input.rmsk} - - {output.map} - 2> {log}"


# SQuIRE clean folder: one *.bed plus one *_copies.txt.
# Count only checks for these two files, so `squire Build` (which would
# download the current UCSC rmsk and merge fragments with bedtools) is skipped.
rule squire_clean_folder:
    input:
        map = f"{OUTDIR}/annotation/te_locus_map_full.tsv",
    output:
        bed    = f"{OUTDIR}/annotation/squire_clean/hg38_all.bed",
        copies = f"{OUTDIR}/annotation/squire_clean/hg38_all_copies.txt",
    log:
        f"{OUTDIR}/annotation/squire_clean_folder.log",
    params:
        outdir = f"{OUTDIR}/annotation/squire_clean",
        script = "workflow/scripts/make_squire_clean.py",
    shell:
        "mkdir -p {params.outdir} && "
        "python {params.script} {input.map} {params.outdir} 2> {log}"


# Cross-annotation audit: quantify once, write it down, never argue again.
rule annotation_check:
    input:
        locus_map  = f"{OUTDIR}/annotation/te_locus_map.tsv",
        squire_bed = f"{OUTDIR}/annotation/squire_clean/hg38_all.bed",
        provirus   = P["provirus_gtf"],
    output:
        f"{OUTDIR}/annotation/annotation_check.tsv",
    params:
        script = "workflow/scripts/annotation_check.py",
    log:
        f"{OUTDIR}/annotation/annotation_check.log",
    shell:
        "mkdir -p $(dirname {output}) && "
        "python {params.script} {input.locus_map} {input.squire_bed} "
        "{input.provirus} {output} 2> {log}"
