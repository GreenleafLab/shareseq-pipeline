#!/bin/bash
set -euo pipefail

# Entry point for processing already-anonymized SHARE-seq FASTQs.
# Runs ingest_anonymized.smk (concatenate SRA parts, count reads) then shareseq.smk.

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
    snakemake --profile="$(pwd)/profile" -s ingest_anonymized.smk --configfile "$configfile" --notemp --use-singularity "${extra_args[@]}"
    snakemake --profile="$(pwd)/profile" -s shareseq.smk          --configfile "$configfile" --notemp --use-singularity "${extra_args[@]}"
else
    echo "Running in local mode"
    snakemake --profile="$(pwd)/profile" -s ingest_anonymized.smk --configfile "$configfile" --notemp "${extra_args[@]}"
    snakemake --profile="$(pwd)/profile" -s shareseq.smk          --configfile "$configfile" --notemp "${extra_args[@]}"
fi
