#!/bin/bash
set -euo pipefail

# If container given in config, run in container mode
# Keep all temporary files for now
if [ ! -z "$(grep -v ^\# $1 | grep -e 'singularity:' | cut -d' ' -f2)" ]
then
    container=$(grep -v ^\# $1 | grep -e 'singularity:' | cut -d' ' -f2)
    echo "Running in singularity container: $container"
    snakemake --profile=$(pwd)/profile -s prep_fastq.smk --configfile $1 --notemp --use-singularity
    snakemake --profile=$(pwd)/profile -s shareseq.smk --configfile $1 --notemp --use-singularity 
else
    echo "Running in local mode"
    snakemake --profile=$(pwd)/profile -s prep_fastq.smk --configfile $1 --notemp
    snakemake --profile=$(pwd)/profile -s shareseq.smk --configfile $1 --notemp
fi
