#!/bin/bash
set -euo pipefail
smk_output=$(grep 'output_dir:' $1 | cut -d: -f2 | cut -d\" -f2)

# overview
Rscript scripts/plotting/plot_overview.R $smk_output plots/overview/

# ATAC
if [ -d $smk_output/ATAC ] ; then
    Rscript scripts/plotting/plot_atac.R $smk_output plots/atac/ hg38
fi

# RNA
if [ -d $smk_output/RNA ] ; then
    Rscript scripts/plotting/plot_rna.R $smk_output plots/rna/ F
fi

# combine outputs into a single pdf
allpdf=$(ls plots/*/*.pdf)
pdfunite $allpdf plots/summary.pdf

