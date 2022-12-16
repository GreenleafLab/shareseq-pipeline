#!/bin/bash
set -euo pipefail

# Keep all temporary files for now
snakemake --profile=$(pwd)/profile -s prep_fastq.smk --configfile $1 --notemp

snakemake --profile=$(pwd)/profile -s shareseq.smk --configfile $1 --notemp
