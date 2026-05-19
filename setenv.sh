#!/bin/bash
# Environment setup for shareseq-pipeline runs on the process-anonymize branch.
# Source (don't execute) before invoking snakemake on this project:
#   source setenv.sh
#
conda activate shareseq

ml biology bcl2fastq/2.20 R/4.1.2 samtools/1.16.1 sra-tools/3.0.7
ml star/2.5.4b
ml star/2.7.10b   # must be loaded after 2.5.4b so it wins on PATH

ml physics geos
ml poppler

ml java/17.0.4    # for variant bams
