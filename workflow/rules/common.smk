# -----------------------------------------------------------------------------
# Shared variables and functions
# -----------------------------------------------------------------------------
import os
import math
import zlib
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

# -----------------------------------------------------------------------------
# Temporary-file placement
#
# Each of the three tools writes GB-scale scratch files. Putting all of them on
# one HDD at high concurrency collapses that disk to ~46 MB/s (measured
# 2026-10-28). config["tmp_roots"] lists candidate roots with a share and the
# Snakemake resource name that limits how many jobs may use that root at once.
# The assignment is a hash of the sample name, so it does not change between
# runs and a resumed run puts a sample in the same place.
# -----------------------------------------------------------------------------
TMP_ROOTS = config.get("tmp_roots", [])
TMP_GB = config.get("resources", {}).get("tmp_gb_per_job", 6)


def tmp_root_index(sample):
    if not TMP_ROOTS:
        return None
    weights = [max(1, int(round(r["share"] * 1000))) for r in TMP_ROOTS]
    h = zlib.crc32(sample.encode()) % sum(weights)
    acc = 0
    for i, w in enumerate(weights):
        acc += w
        if h < acc:
            return i
    return len(TMP_ROOTS) - 1


def tmp_dir(wildcards):
    idx = tmp_root_index(wildcards.sample)
    if idx is None:
        return f"{OUTDIR}/_tmp/{wildcards.sample}"
    return f"{TMP_ROOTS[idx]['path']}/{wildcards.sample}"


def tmp_gb(wildcards):
    """Charge this job's temp budget to its root's resource counter only.

    A flat budget is wrong: measured usage ranges from 3.5 GB (small BAM) to
    16.6 GB (5 GB BAM, telescope collate + TEcount sort alive together), so a
    flat 6 GB per job over-subscribed /dev/shm and the tmpfs filled up, which
    made intermediate bed files disappear mid-run. Temp size tracks the BAM
    size, so it is estimated from it: ~3.5x the input BAM, floor 6 GB.
    """
    idx = tmp_root_index(wildcards.sample)
    if idx is None:
        return {}
    return {TMP_ROOTS[idx]["resource"]: sample_tmp_gb(wildcards.sample)}


_BAM_GB = {}


def sample_tmp_gb(sample):
    if sample not in _BAM_GB:
        bam = next((s["bam"] for s in SAMPLES if s["sample"] == sample), None)
        gb = os.path.getsize(bam) / 2**30 if bam and os.path.exists(bam) else 2.0
        _BAM_GB[sample] = max(6, int(math.ceil(3.5 * gb)))
    return _BAM_GB[sample]


def tmp_shm_gb(wildcards):
    return tmp_gb(wildcards).get("shm_gb", 0)


def tmp_hdd_gb(wildcards):
    return tmp_gb(wildcards).get("hdd_gb", 0)


def tmp_hdd2_gb(wildcards):
    return tmp_gb(wildcards).get("hdd2_gb", 0)

def gene_count_outputs():
    """Per-sample featureCounts outputs, used as the merge rule's input."""
    if not config["tools"].get("gene_counts", {}).get("run"):
        return []
    return [f"{OUTDIR}/gene_counts/{s['cohort']}/{s['sample']}/"
            f"{s['sample']}_featureCounts.txt" for s in SAMPLES]


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
        # ERVmap is imported, not rerun, but its normalised table must exist
        # before aggregation, otherwise the ERVmap column silently stays empty.
        # tools.ervmap.run=false means "do not re-quantify ERVmap"; importing the
        # existing raw results is still required.
        raw = f"{OUTDIR}/ervmap_raw/{c}/{smp}ERVresults.txt"
        if os.path.exists(raw):
            out.append(f"{OUTDIR}/ervmap/{c}/{smp}/{smp}.locus.tsv")
    return out
