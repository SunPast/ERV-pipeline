# -----------------------------------------------------------------------------
# FASTQ -> BAM
#
#   fastp -> STAR -> samtools sort/index
#
# A sample is aligned here when its `fq1` column is filled in samples.tsv;
# otherwise the `bam` column is used as-is (legacy cohorts).
# Both paired-end (fq1+fq2) and single-end (fq1 only) are supported.
#
# Outputs
#   fastp/{cohort}/{sample}/{sample}.trim.R1.fq.gz   (+ .trim.R2.fq.gz when paired)
#   star/{cohort}/{sample}/Aligned.out.bam, {sample}_Log.final.out
#   bam/{cohort}/{sample}/{sample}.sorted.bam(+.bai)   <- consumed by every tool rule
# -----------------------------------------------------------------------------
_FASTP = config.get("tools", {}).get("fastp", {})
_STAR = config.get("tools", {}).get("star", {})
_FASTP_RUN = _FASTP.get("run", True)


def _sample(wildcards):
    return next(s for s in SAMPLES if s["sample"] == wildcards.sample)


def is_paired(wildcards):
    return bool(_sample(wildcards).get("fq2"))


def _fastp_dir(wc):
    return f"{OUTDIR}/fastp/{wc.cohort}/{wc.sample}"


def trimmed_r1(wildcards):
    return f"{_fastp_dir(wildcards)}/{wildcards.sample}.trim.R1.fq.gz"


def trimmed_r2(wildcards):
    return f"{_fastp_dir(wildcards)}/{wildcards.sample}.trim.R2.fq.gz"


def star_inputs(wildcards):
    """Read files for STAR; keeps the fastp rule in the DAG."""
    s = _sample(wildcards)
    if _FASTP_RUN:
        return [trimmed_r1(wildcards)] + ([trimmed_r2(wildcards)] if s["fq2"] else [])
    return [s["fq1"]] + ([s["fq2"]] if s["fq2"] else [])


if _FASTP_RUN:
    rule fastp:
        input:
            reads = lambda wc: (
                [_sample(wc)["fq1"], _sample(wc)["fq2"]]
                if _sample(wc)["fq2"] else [_sample(wc)["fq1"]]
            ),
        output:
            r1 = f"{OUTDIR}/fastp/{{cohort}}/{{sample}}/{{sample}}.trim.R1.fq.gz",
            r2 = f"{OUTDIR}/fastp/{{cohort}}/{{sample}}/{{sample}}.trim.R2.fq.gz",
        threads: _FASTP.get("threads", 8)
        resources:
            mem_mb = _FASTP.get("mem_mb", 4000)
        params:
            bin = P.get("fastp_bin", "/usr/bin"),
            extra = _FASTP.get("extra", ""),
            min_length = _FASTP.get("min_length", 36),
            phred = _FASTP.get("qualified_quality_phred", 20),
            unqual_pct = _FASTP.get("unqualified_percent_limit", 40),
            n_base = _FASTP.get("n_base_limit", 5),
            pe = lambda wc: (
                (f"-I {_sample(wc)['fq2']} -O "
                 f"{OUTDIR}/fastp/{wc.cohort}/{wc.sample}/{wc.sample}.trim.R2.fq.gz "
                 "--detect_adapter_for_pe")
                if _sample(wc)["fq2"] else ""
            ),
            html = lambda wc: f"{_fastp_dir(wc)}/{wc.sample}.fastp.html",
            json = lambda wc: f"{_fastp_dir(wc)}/{wc.sample}.fastp.json",
        log:
            f"{OUTDIR}/fastp/{{cohort}}/{{sample}}/fastp.log"
        shell:
            """
            mkdir -p $(dirname {log:q})
            {params.bin}/fastp -w {threads} \
                -i {input.reads[0]:q} -o {output.r1:q} {params.pe} \
                -h {params.html:q} -j {params.json:q} \
                --length_required {params.min_length} \
                --qualified_quality_phred {params.phred} \
                --unqualified_percent_limit {params.unqual_pct} \
                --n_base_limit {params.n_base} \
                {params.extra} > {log} 2>&1
            test -s {output.r2:q} || touch {output.r2:q}
            """


rule star_align:
    input:
        reads = star_inputs,
        index = P["star_index"],
    output:
        bam = f"{OUTDIR}/star/{{cohort}}/{{sample}}/Aligned.out.bam",
        log_final = f"{OUTDIR}/star/{{cohort}}/{{sample}}/{{sample}}_Log.final.out",
    threads: _STAR.get("threads", 12)
    resources:
        mem_mb = _STAR.get("mem_mb", 40000)
    params:
        bin = P.get("star_bin", "/usr/bin"),
        outdir = lambda wc: f"{OUTDIR}/star/{wc.cohort}/{wc.sample}",
        prefix = lambda wc: f"{OUTDIR}/star/{wc.cohort}/{wc.sample}/{wc.sample}_",
        mismatch_max = _STAR.get("out_filter_mismatch_nmax", 15),
        bam_sort_ram = _STAR.get("limit_bam_sort_ram", 39050942993),
        extra = _STAR.get("extra", ""),
    log:
        f"{OUTDIR}/star/{{cohort}}/{{sample}}/star.log"
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log:q})
        {params.bin}/STAR --runThreadN {threads} --genomeDir {input.index} \
            --readFilesIn {input.reads} --readFilesCommand zcat \
            --outSAMunmapped Within \
            --outFilterMismatchNmax {params.mismatch_max} \
            --outSAMattributes All \
            --limitBAMsortRAM {params.bam_sort_ram} \
            --outSAMtype BAM Unsorted \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMstrandField intronMotif \
            --twopassMode Basic \
            --quantMode GeneCounts \
            --outFileNamePrefix {params.prefix} \
            {params.extra} >> {log} 2>&1
        mv -f {params.prefix}Aligned.out.bam {output.bam:q}
        """


rule samtools_sort:
    input:
        bam = f"{OUTDIR}/star/{{cohort}}/{{sample}}/Aligned.out.bam",
    output:
        bam = f"{OUTDIR}/bam/{{cohort}}/{{sample}}/{{sample}}.sorted.bam",
        bai = f"{OUTDIR}/bam/{{cohort}}/{{sample}}/{{sample}}.sorted.bam.bai",
    threads: 8
    resources:
        mem_mb = 8000
    params:
        bin = P.get("samtools_bin", "/usr/bin"),
        tmp = lambda wc: f"{tmp_dir(wc)}_sort",
    log:
        f"{OUTDIR}/bam/{{cohort}}/{{sample}}/sort.log"
    shell:
        """
        mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
        {params.bin}/samtools sort -@ {threads} -m 2G -T {params.tmp} \
            -o {output.bam:q} {input.bam:q} >> {log} 2>&1
        {params.bin}/samtools index -@ 4 {output.bam:q}
        rm -rf {params.tmp}.*
        """
