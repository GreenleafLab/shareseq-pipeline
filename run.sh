#!/bin/bash
ml biology bcl2fastq samtools
set -euo pipefail

# First argument: config file
configfile="$1"
shift  # remove the first argument, leaving all extras in "$@"

# Optional extra args (e.g. -n, --cores, --rerun-incomplete, etc.)
extra_args=("$@")

# If container given in config, run in container mode
# Keep all temporary files for now
if [ ! -z "$(grep -v ^\# $configfile | grep -e 'singularity:' | cut -d' ' -f2)" ]
then
    container=$(grep -v ^\# $configfile | grep -e 'singularity:' | cut -d' ' -f2)
    echo "Running in singularity container: $container"
    snakemake --profile="$(pwd)/profile" -s prep_fastq.smk --configfile "$configfile" --notemp --use-singularity "${extra_args[@]}"
    snakemake --profile="$(pwd)/profile" -s anonymize.smk --configfile "$configfile" --notemp --use-singularity "${extra_args[@]}"
else
    echo "Running in local mode"
    snakemake --profile="$(pwd)/profile" -s prep_fastq.smk --configfile "$configfile" --notemp "${extra_args[@]}"
    snakemake --profile="$(pwd)/profile" -s anonymize.smk --configfile "$configfile" --notemp "${extra_args[@]}"
fi
