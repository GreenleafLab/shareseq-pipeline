# Ingestion of already-anonymized SHARE-seq FASTQs (e.g. from SRA) for shareseq.smk
# Author: Selin Jessa
# Last Modified: 05/14/2026

# Replaces prep_fastq.smk for inputs that are already demultiplexed per sample
# and that already carry the SHARE-seq index info `1:N:0:[I1]+[I2]` in the read
# header (produced by anonymize.smk -> merge_samples and redistributed via SRA).
# 1. Split each per-sample merged FASTQ by I2 into per-sublibrary `.fastq.zst`
#    files (mirroring the per-sublibrary structure prep_fastq.smk produces for
#    bcl inputs). The I2 -> sublibrary map is config["sequencing"][run]["sublibraries"][assay].
# 2. Count the number of reads per staged fastq file (per-sublibrary granularity).

import json
import os
import sys

import utils

workdir: config["output_dir"]

# global singularity container to use
# only set if container given in config and is not none
if "singularity" in config.keys() and config["singularity"]:
    singularity: config["singularity"]

wildcard_constraints:
    sequencing_path = "(ATAC|RNA)/[^/]+/[^/]+",
    assay = "ATAC|RNA",
    sample = "[^/]+",
    read = "R[12]",

utils.string_only_keys(config)

def get_anon_input_fastq(assay, sample, read, config):
    """Return the single merged per-sample anon FASTQ path produced by prep_sra.smk
    (or already pre-staged by the user)."""
    for run_id, run in config["sequencing"].items():
        if run.get("type") != "anonymized":
            continue
        if sample not in (run.get(f"{assay}_samples") or []):
            continue
        data_dir = run["data_dir"].rstrip("/")
        path = os.path.join(data_dir, f"{assay}_{sample}_anon_{read}.fastq.gz")
        if not os.path.exists(path):
            raise RuntimeError(
                f"Missing anon FASTQ: {path}. Run prep_sra.smk first or pre-stage the file."
            )
        return path
    raise RuntimeError(f"No anonymized sequencing run lists {sample} under {assay}_samples")


def get_single_anon_sample(run_id, assay):
    """For an anonymized run + assay, return the single SampleID listed there.
    The split-by-I2 rule below assumes exactly one sample per (run, assay)."""
    samples = config["sequencing"][run_id].get(f"{assay}_samples") or []
    if len(samples) != 1:
        raise RuntimeError(
            f"split_anon expects exactly one sample per (run, assay); "
            f"{run_id}/{assay} has: {samples}"
        )
    return samples[0]


anon_paths = (
    utils.get_sequencing_paths("ATAC", config, run_types=["anonymized"]) +
    utils.get_sequencing_paths("RNA", config, run_types=["anonymized"])
)

outputs = (
    expand("{sequencing_path}/read_count.txt", sequencing_path=anon_paths) +
    expand("staged_fastq/{sequencing_path}_{read}.fastq.zst", sequencing_path=anon_paths, read=["R1", "R2"])
)

if "filter_dag" in config.keys() and config["filter_dag"] == "false":
    filtered_outputs = outputs
else:
    filtered_outputs = []
    for o in outputs:
        if os.path.exists(o):
            print(f"Skipping existing output: {o}", file=sys.stderr)
        else:
            filtered_outputs.append(o)

localrules: all
rule all:
    input: filtered_outputs


#############################
### Split per-sample merged anon FASTQ by I2 into per-sublibrary files.
### One rule per (assay, run_id, sample) -- generated in a Python loop so each
### rule statically declares its (config-driven) sublibrary outputs.
#############################

def _split_rule_name(assay, run_id, sample):
    return f"split_anon_{assay}_{run_id}_{sample}".replace("-", "_").replace(".", "_")


for _run_id, _run in config["sequencing"].items():
    if _run.get("type") != "anonymized":
        continue
    _data_dir = _run["data_dir"].rstrip("/")
    for _assay in ("ATAC", "RNA"):
        _sublibs = (_run.get("sublibraries") or {}).get(_assay) or {}
        _samples = _run.get(f"{_assay}_samples") or []
        if not (_sublibs and _samples):
            continue
        for _sample in _samples:
            _sublib_names = list(_sublibs.keys())
            _outputs = [
                f"staged_fastq/{_assay}/{_run_id}/{_sublib}_{_read}.fastq.zst"
                for _sublib in _sublib_names for _read in ("R1", "R2")
            ]
            _stats_path = f"staged_fastq/logs/{_assay}_{_run_id}_{_sample}_split_stats.json"
            _log_path = f"staged_fastq/logs/{_assay}_{_run_id}_{_sample}_split.log"
            _i2_map_json = json.dumps(_sublibs)
            _outdir = f"staged_fastq/{_assay}/{_run_id}"
            _r1_in = f"{_data_dir}/{_assay}_{_sample}_anon_R1.fastq.gz"
            _r2_in = f"{_data_dir}/{_assay}_{_sample}_anon_R2.fastq.gz"

            rule:
                name: _split_rule_name(_assay, _run_id, _sample)
                input:
                    R1 = _r1_in,
                    R2 = _r2_in,
                output:
                    fastqs = _outputs,
                    stats = _stats_path,
                log: _log_path
                params:
                    script = "scripts/shareseq/split_anon_fastq_by_i2.py",
                    i2_map = _i2_map_json,
                    outdir = _outdir,
                resources:
                    runtime = 6 * 60,
                    mem_mb = 16_000,
                threads: 8
                shell:
                    "python {params.script} "
                    " --r1-in {input.R1} --r2-in {input.R2} "
                    " --i2-map '{params.i2_map}' "
                    " --out-dir {params.outdir} "
                    " --stats {output.stats} "
                    " --threads {threads} "
                    " 2> {log}"


#############################
### Count reads (per-sublibrary)
#############################

rule count_reads:
    input: lambda w: utils.fastq_path(w.sequencing_path, "R1", config)
    output: "{sequencing_path}/read_count.txt"
    params:
        decompress = lambda w: utils.fastq_decompress(w.sequencing_path, config),
        truncate_test_chunks = lambda w: f" | head -n {config['chunk_size']*config['test_chunks']*4} " if "test_chunks" in config else ""
    resources:
        runtime = 5 * 60,
    shell: "{params.decompress} {input} {params.truncate_test_chunks} | awk -c 'END{{print int(NR/4)}}' > {output}"
           "|| if [[ $? -eq 141 ]]; then true; else exit $?; fi" # Ignore spurious exit code 141 produced with test_chunks set
