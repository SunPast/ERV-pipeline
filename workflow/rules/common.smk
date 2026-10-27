# -----------------------------------------------------------------------------
# Shared variables and functions
# -----------------------------------------------------------------------------
import os
import pandas as pd

OUTDIR = config["outdir"]
P = config["paths"]

SAMPLES_FILE = config.get("samples", "config/samples.example.tsv")
_samples = pd.read_csv(SAMPLES_FILE, sep="\t", dtype=str).fillna("")
SAMPLES = _samples.to_dict("records")

def sdir(tool, cohort, sample):
    return f"{OUTDIR}/{tool}/{cohort}/{sample}"

def cohort_dir(tool, cohort):
    return f"{OUTDIR}/{tool}/{cohort}"

# Sampler: Snakemake 8 uses input functions; here we keep a 7/8-compatible lambda form
def bam_of(wildcards):
    return next(s["bam"] for s in SAMPLES if s["sample"] == wildcards.sample)

def fq_of(wildcards):
    s = next(s for s in SAMPLES if s["sample"] == wildcards.sample)
    return s["fq1"], s["fq2"]

def run_env():
    """Export the bin directory containing telescope/samtools (same env)"""
    return f"export PATH={P['telescope_bin']}:$PATH && "

def apptainer():
    return f"{P['apptainer_bin']}/apptainer"

def all_tool_outputs():
    """Per-sample outputs of the three rerun tools, used to wire the aggregate step.

    Only the tools that are rerun here are listed. ERVmap is imported from
    existing results and is therefore not a dependency of `all`.
    """
    out = []
    for s in SAMPLES:
        c, smp = s["cohort"], s["sample"]
        if config["tools"]["tetranscripts"]["run"]:
            out.append(f"{OUTDIR}/tetranscripts/{c}/{smp}/{smp}.cntTable")
        if config["tools"]["telescope"]["run"]:
            out.append(f"{OUTDIR}/telescope/{c}/{smp}/{smp}.report.tsv")
        if config["tools"]["squire"]["run"]:
            out.append(f"{OUTDIR}/squire/{c}/{smp}/{smp}_subFcounts.txt")
    return out
