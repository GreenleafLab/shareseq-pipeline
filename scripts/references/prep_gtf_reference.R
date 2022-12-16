suppressPackageStartupMessages({
    library(stringr)
    library(readr)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

input_path <- args[1] # Can be a URL or a file path to a gtf
output_path <- args[2]

gtf_colnames <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
gtf <- readr::read_tsv(input_path, comment = "#", col_names = gtf_colnames, col_types = "ccciicccc", progress=TRUE)

gene_id <- stringr::str_match(gtf$attributes, 'gene_id "([^"]*)"')[,2]
gene_type <- stringr::str_match(gtf$attributes, 'gene_type "([^"]*)"')[,2]

keeper_types <- "lncRNA|protein_coding|IG_.*|TR_.*"
keeper_genes <- unique(gene_id[stringr::str_detect(gene_type, keeper_types)])

gtf_filt <- gtf[gene_id %in% keeper_genes,]

header_lines <- c(
    sprintf("##description: Filtered gtf from: %s", input_path),
    sprintf("##filtering: gene_type matches: %s", keeper_types)
)

write_lines(header_lines, output_path)
write_tsv(gtf_filt, output_path, append=TRUE, col_names=FALSE, quote="none", escape="none")