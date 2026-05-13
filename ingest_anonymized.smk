# Ingestion of already-anonymized SHARE-seq FASTQs (e.g. from SRA) for shareseq.smk
# Author: Selin Jessa
# Last Modified: 05/13/2026

# Replaces prep_fastq.smk for inputs that are already demultiplexed per sample
# and that already carry the SHARE-seq index info `1:N:0:[I1]+[I2]` in the read
# header (produced by anonymize.smk -> merge_samples and redistributed via SRA).
# 1. Concatenate any SRA-chunked .part_*.fastq.gz files for each sample/read into
#    a single .fastq.zst, mirroring the prep_fastq.smk bcl2fastq output format
# 2. Count the number of reads per staged fastq file

import os
import glob
import sys

import utils

workdir: config["output_dir"]

# global singularity container to use
# only set if container given in config and is not none
if "singularity" in config.keys() and config["singularity"]:
    singularity: config["singularity"]

wildcard_constraints:
    sequencing_path = "(ATAC|RNA)/samples/[^/]+",
    assay = "ATAC|RNA",
    sample = "[^/]+",
    read = "R[12]",

utils.string_only_keys(config)

def get_anon_parts(assay, sample, read, config):
    """Return the sorted list of source FASTQ paths in data_dir for a sample+read.

    Accepts either a single file `{ASSAY}_{Sample}_anon_{read}.fastq.gz` or a
    set of SRA-chunked parts `{ASSAY}_{Sample}_anon_{read}.part_*.fastq.gz`.
    """
    for run_id, run in config["sequencing"].items():
        if run.get("type") != "anonymized":
            continue
        if sample not in (run.get(f"{assay}_samples") or []):
            continue
        data_dir = run["data_dir"]
        single = os.path.join(data_dir, f"{assay}_{sample}_anon_{read}.fastq.gz")
        if os.path.exists(single):
            return [single]
        parts = sorted(glob.glob(os.path.join(data_dir, f"{assay}_{sample}_anon_{read}.part_*.fastq.gz")))
        if parts:
            return parts
        raise RuntimeError(
            f"No FASTQ files found in {data_dir} for {assay}/{sample}/{read} "
            f"(looked for '{assay}_{sample}_anon_{read}.fastq.gz' and "
            f"'{assay}_{sample}_anon_{read}.part_*.fastq.gz')"
        )
    raise RuntimeError(f"No anonymized sequencing run lists {sample} under {assay}_samples")

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
### Concatenate SRA-chunked parts into a single .fastq.zst per sample/read
#############################

rule concatenate_anon_fastqs:
    input:
        fastqs = lambda w: get_anon_parts(w.assay, w.sample, w.read, config)
    output:
        fastq = "staged_fastq/{assay}/samples/{sample}_{read}.fastq.zst"
    resources:
        runtime = 5 * 60,
    threads: 4
    log: "staged_fastq/logs/{assay}_{sample}_{read}_concat.log"
    # SRA-anonymized FASTQs from anonymize.smk are written by `samtools fastq`,
    # which separates the read ID from the `1:N:0:...` index info with a TAB.
    # shareseq.smk's match_barcodes.py splits on a single space (matches
    # bcl2fastq output); the awk pass normalizes the TAB on header lines so the
    # downstream pipeline behaves identically to a fresh bcl2fastq run.
    shell: "zcat {input.fastqs} | awk 'NR%4==1 {{sub(\"\\t\", \" \")}} 1' | zstd --fast -q -T{threads} -o {output.fastq} 2> {log}"

#############################
### Count reads
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
