# SRA fetch + extract stage for the anonymized-ingest workflow
# Author: Selin Jessa
# Last Modified: 05/14/2026
#
# Optional pre-stage for ingest_anonymized.smk. For samples whose anonymized
# FASTQs live on SRA, declare a {sample: SRR} map under
# `sequencing.{run_id}.{ATAC,RNA}_sra` in the config; this snakefile fetches
# each SRR with prefetch + fasterq-dump (multi-threaded, ~5x faster than
# fastq-dump), reconstructs the original Illumina header from SPOT_GROUP via
# --seq-defline, then pigz-compresses the output and lands it in `data_dir`
# under the naming scheme ingest_anonymized.smk expects
# ({ASSAY}_{Sample}_anon_{R1|R2}.fastq.gz).
#
# Samples not listed in {ATAC,RNA}_sra are silently skipped (assumed to be
# pre-staged locally). If no run declares ATAC_sra/RNA_sra at all, this
# snakefile is a no-op and behaves identically to run_anon.sh.
#
# Requires sra-tools (prefetch, fasterq-dump) + pigz on $PATH. The user is
# expected to have sourced an env script (e.g. the project-level setenv.sh
# that does `conda activate shareseq` + `ml biology sra-tools/3.0.7`) before
# invoking snakemake. Same convention as prep_fastq.smk / anonymize.smk.
#
# Disk usage: fasterq-dump spills temp files of ~10x the .sra size during
# extraction, so a 35 GB scATAC SRA needs ~350 GB of temp space. We point
# --temp at $L_SCRATCH (Sherlock local-NVMe scratch); fall back to /tmp if
# unset. The final gzipped output is roughly 1.5x the .sra size.

import os
import sys

import utils

workdir: config["output_dir"]

# global singularity container to use
# only set if container given in config and is not none
if "singularity" in config.keys() and config["singularity"]:
    singularity: config["singularity"]

wildcard_constraints:
    data_dir = ".+",
    assay = "ATAC|RNA",
    sample = "[^/]+",

utils.string_only_keys(config)


def iter_sra_targets(config):
    """Yield (run_id, assay, sample, srr, data_dir) tuples for SRA-declared samples."""
    for run_id, run in config["sequencing"].items():
        if run.get("type") != "anonymized":
            continue
        data_dir = run.get("data_dir")
        if not data_dir:
            continue
        for assay in ("ATAC", "RNA"):
            sra_map = run.get(f"{assay}_sra") or {}
            for sample, srr in sra_map.items():
                yield run_id, assay, sample, srr, data_dir


def lookup_srr(wildcards):
    """Resolve the SRR for a (data_dir, assay, sample) wildcard match."""
    target_dir = wildcards.data_dir.rstrip("/")
    matches = [
        srr for _, assay, sample, srr, data_dir in iter_sra_targets(config)
        if assay == wildcards.assay
        and sample == wildcards.sample
        and data_dir.rstrip("/") == target_dir
    ]
    if not matches:
        raise RuntimeError(
            f"No SRR mapping found for {wildcards.assay}/{wildcards.sample} in "
            f"data_dir '{wildcards.data_dir}'. Check sequencing.*.{wildcards.assay}_sra."
        )
    if len(set(matches)) > 1:
        raise RuntimeError(
            f"Multiple distinct SRR mappings for {wildcards.assay}/{wildcards.sample} "
            f"in data_dir '{wildcards.data_dir}': {sorted(set(matches))}"
        )
    return matches[0]


outputs = [
    f"{data_dir.rstrip('/')}/{assay}_{sample}_anon_{read}.fastq.gz"
    for (_, assay, sample, _, data_dir) in iter_sra_targets(config)
    for read in ("R1", "R2")
]

if "filter_dag" in config.keys() and config["filter_dag"] == "false":
    filtered_outputs = outputs
else:
    filtered_outputs = []
    for o in outputs:
        if os.path.exists(o):
            print(f"Skipping existing SRA output: {o}", file=sys.stderr)
        else:
            filtered_outputs.append(o)

localrules: all
rule all:
    input: filtered_outputs


#############################
### Fetch and extract one SRR -> {ASSAY}_{sample}_anon_{R1,R2}.fastq.gz
#############################

rule fetch_and_extract_sra:
    output:
        R1 = "{data_dir}/{assay}_{sample}_anon_R1.fastq.gz",
        R2 = "{data_dir}/{assay}_{sample}_anon_R2.fastq.gz",
    params:
        srr = lookup_srr,
        workdir = lambda w: f"{w.data_dir.rstrip('/')}/_sra_tmp/{w.assay}_{w.sample}",
    resources:
        runtime = 4 * 60,    # fasterq-dump is ~5x faster than fastq-dump
    threads: 8
    log: "{data_dir}/logs/{assay}_{sample}_sra.log"
    shell:
        """
        mkdir -p {params.workdir}
        tempdir="${{L_SCRATCH:-/tmp}}/fasterq-dump_{params.srr}_$$"
        mkdir -p "$tempdir"
        trap 'rm -rf "$tempdir"' EXIT
        ( cd {params.workdir} && \
            prefetch --max-size 100g {params.srr} && \
            fasterq-dump --split-files \
                --seq-defline '@$sn 1:N:0:$sg' --qual-defline '+' \
                -e {threads} -t "$tempdir" \
                {params.srr} && \
            pigz -p {threads} {params.srr}_1.fastq {params.srr}_2.fastq ) 2> {log}
        mv {params.workdir}/{params.srr}_1.fastq.gz {output.R1}
        mv {params.workdir}/{params.srr}_2.fastq.gz {output.R2}
        rm -rf {params.workdir}
        """
