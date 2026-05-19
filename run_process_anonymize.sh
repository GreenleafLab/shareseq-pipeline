#!/bin/bash
set -euo pipefail

# Entry point for processing already-anonymized SHARE-seq FASTQs, with an
# optional SRA-download pre-stage. Runs prep_sra.smk (fetch any SRRs declared
# via ATAC_sra/RNA_sra into data_dir), then ingest_anonymized.smk (concatenate
# parts, count reads), then shareseq.smk.
#
# prep_sra.smk is a no-op when no run in the config declares ATAC_sra/RNA_sra,
# so this wrapper is also a drop-in replacement for run_anon.sh for runs whose
# fastqs are already local.
#
# Expects the shareseq env to be active before invocation (e.g. source the
# project's setenv.sh, which does conda activate + ml loads, including
# sra-tools/3.0.7 for the prep_sra step). Same expectation as run_anon.sh.

# First argument: config file
configfile="$1"
shift || true  # remove the first argument, safely even if none follow

# Optional extra args (e.g. -n, --cores, etc.)
extra_args=("$@")
# If none given, default to an empty array
if [ ${#extra_args[@]} -eq 0 ]; then
    extra_args=("")
fi

# If container given in config, run in container mode
# Keep all temporary files for now
if [ ! -z "$(grep -v ^\# $configfile | grep -e 'singularity:' | cut -d' ' -f2)" ]
then
    container=$(grep -v ^\# $configfile | grep -e 'singularity:' | cut -d' ' -f2)
    echo "Running in singularity container: $container"
    snakemake --profile="$(pwd)/profile" -s prep_sra.smk          --configfile "$configfile" --notemp --use-singularity "${extra_args[@]}"
    snakemake --profile="$(pwd)/profile" -s ingest_anonymized.smk --configfile "$configfile" --notemp --use-singularity "${extra_args[@]}"
    snakemake --profile="$(pwd)/profile" -s shareseq.smk          --configfile "$configfile" --notemp --use-singularity "${extra_args[@]}"
else
    echo "Running in local mode"
    snakemake --profile="$(pwd)/profile" -s prep_sra.smk          --configfile "$configfile" --notemp "${extra_args[@]}"
    snakemake --profile="$(pwd)/profile" -s ingest_anonymized.smk --configfile "$configfile" --notemp "${extra_args[@]}"
    snakemake --profile="$(pwd)/profile" -s shareseq.smk          --configfile "$configfile" --notemp "${extra_args[@]}"
fi
