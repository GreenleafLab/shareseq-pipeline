suppressPackageStartupMessages({
  library(BPCells)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
})


########## arguments ############
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)

input_path <- args[1] # path to the shareseq.smk output folder
output_path <- args[2] # path to save script outputs (plots, tables, objects)
genome <- args[3] # only hg38 support for now

stopifnot(dir.exists(paste0(input_path, "/ATAC")))
stopifnot(genome=="hg38") # TODO support for genomes other than hg38

########## functions ############
get_references <- function(output_path, genome){
  # Reference annotations
  refdir <- paste0(output_path, "/references")
  genes <- read_gencode_transcripts( 
    refdir, 
    release="42", 
    transcript_choice="MANE_Select",
    annotation_set = "basic", 
    features="transcript" # Make sure to set this so we don't get exons as well
  )
  
  blacklist <- read_encode_blacklist(refdir, genome=genome)
  chrom_sizes <- read_ucsc_chrom_sizes(refdir, genome=genome)
  return(list(genes=genes, blacklist=blacklist, chrom_sizes=chrom_sizes))
}

 
plot_frag_tss_persublib <- function(input_path, output_path, genes, blacklist){
  tmp <- list.files(paste0(input_path,"/ATAC/sublibraries"), recursive=T)
  frag.file.list <- tmp[grep("\\/fragments.tsv.gz.tbi", tmp)] %>% 
    strsplit(".tbi") %>% unlist
  
  pdf(paste0(output_path,"/ATAC_fragment_TSS_profile_persublib.pdf"), width=10, height=4)
  for (frag.file in frag.file.list) {
    sample <- gsub("/", "_",  frag.file %>% strsplit("/fragments.tsv.gz") %>% unlist)
    fragdir <- paste0(output_path, "/frags_", sample)
    fragfile <- paste0(input_path, "/ATAC/sublibraries/", frag.file)
    
    message(fragfile)
    message("reading data")
    # Check if we already ran import
    if (!file.exists(fragdir)) {
      frags_raw <- open_fragments_10x(fragfile) %>%
        write_fragments_dir(fragdir)
    } else {
      frags_raw <- open_fragments_dir(fragdir)
    }
    
    # Calculate ATAC-seq quality-control metrics
    message("calculating ATAC qc")
    atac_qc <- qc_scATAC(frags_raw, genes, blacklist)
    
    message("plotting qc plots")
    p <- plot_fragment_length(frags_raw) + ggtitle(sample) + plot_tss_profile(frags_raw, genes) + plot_tss_scatter(atac_qc, min_frags=1000, min_tss=10)
    print(p)
  }
  invisible(dev.off())
}

plot_frag_tss_persample <- function(input_path, output_path, genes, blacklist){
  tmp <- list.files(paste0(input_path,"/ATAC/samples"), recursive=F)
  frag.file.list <- tmp[grep(".fragments.tsv.gz.tbi", tmp)] %>% 
    strsplit(".tbi") %>% unlist
  
  pdf(paste0(output_path,"/ATAC_fragment_TSS_profile_persample.pdf"), width=10, height=4)
  for (frag.file in frag.file.list) {
    sample <- gsub("/", "_",  frag.file %>% strsplit(".fragments.tsv.gz") %>% unlist)
    fragdir <- paste0(output_path, "/frags_", sample)
    fragfile <- paste0(input_path, "/ATAC/samples/", frag.file)
    
    message(fragfile)
    message("reading data")
    # Check if we already ran import
    if (!file.exists(fragdir)) {
      frags_raw <- open_fragments_10x(fragfile) %>%
        write_fragments_dir(fragdir)
    } else {
      frags_raw <- open_fragments_dir(fragdir)
    }
    
    # Calculate ATAC-seq quality-control metrics
    message("calculating ATAC qc")
    atac_qc <- qc_scATAC(frags_raw, genes, blacklist)
    
    message("plotting qc plots")
    p <- plot_fragment_length(frags_raw) + ggtitle(sample) + plot_tss_profile(frags_raw, genes) + plot_tss_scatter(atac_qc, min_frags=1000, min_tss=10)
    print(p)
  }
  invisible(dev.off())
}

########## main ############
if (!dir.exists(output_path)){
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
}

refs <- get_references(output_path, genome)
plot_frag_tss_persublib(input_path, output_path, refs$genes, refs$blacklist)
#plot_frag_tss_persample(input_path, output_path, refs$genes, refs$blacklist) # optional
