suppressPackageStartupMessages({
  library(rjson)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})


########## arguments ############
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

input_path <- args[1] # path to the shareseq.smk output folder
output_path <- args[2] # path to save script outputs (plots, tables, objects)


########## functions ############
plot_barcode_stats <- function(bc.stats.data, output_path){
  status.list <- c("exact match", "1bp mismatch in 1 barcode", "1bp mismatch in 2 barcodes", "1bp mismatch in 3 barcodes", "unmatched")
  
  df <- unlist(bc.stats.data['total_mismatch_histogram',]) %>% rbind %>% t %>% as.data.frame
  colnames(df) <- "reads"
  rownames(df) <- NULL
  df[,"matching_status"] <- rep(status.list, dim(bc.stats.data)[2]) %>% factor(levels = rev(status.list))
  df[,"sublibrary"] <- lapply(colnames(bc.stats.data), function(n){rep(n,5)}) %>% unlist
  
  # plot percentage of matched barcodes
  p <- ggplot(df, aes(fill=matching_status, y=reads, x=sublibrary)) + geom_bar(position="fill", stat= "identity") + 
      theme_classic() + theme(text=element_text(family="sans", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      scale_fill_brewer(palette="Blues") 
  ggsave(plot=p, paste0(output_path,"/barcode_stats.pdf"), dpi=600, width=10, height=5)
}


plot_barcode_distr <- function(bc.stats.data, output_path){
  get_bc_counts <- function(bc, bc.stats.data){
    dfout <- data.frame()
    for (sample in colnames(bc.stats.data)){
      df <- unlist(bc.stats.data[bc,sample][[1]]$exact_match) %>% rbind %>% t %>%  as.data.frame
      dff <- data.frame(sample = df[sort(rownames(df)),], row.names=sort(rownames(df)))
      colnames(dff) <- sample
      if (dim(dfout)[2] == 0){
        dfout <- dff
      } else {
        dfout <- cbind(dfout, dff)
      }
    }
    
    dfout$col <- (0:95)%%12 + 1
    dfout$row <- (0:95)%/%12 + 1
    return(dfout)
  }
  
  df1 <- get_bc_counts("BC1", bc.stats.data)
  df2 <- get_bc_counts("BC2", bc.stats.data)
  df3 <- get_bc_counts("BC3", bc.stats.data)
  
  pdf(paste0(output_path,"/barcode_distribution.pdf"), width=10, height=8)
  for (sample in colnames(bc.stats.data)){
    p1 <- ggplot(df1,aes_string(x="col",y="row",col=sample)) + geom_point(aes_string(size=sample))+ theme_classic() +
      scale_color_gradient(low="blue",high="red") + ggtitle("Barcode 1 distribution across plate")
    p2 <- ggplot(df2,aes_string(x="col",y="row",col=sample)) + geom_point(aes_string(size=sample))+ theme_classic() +
      scale_color_gradient(low="blue",high="red") + ggtitle("Barcode 2 distribution across plate")
    p3 <- ggplot(df3,aes_string(x="col",y="row",col=sample)) + geom_point(aes_string(size=sample))+ theme_classic() +
      scale_color_gradient(low="blue",high="red") + ggtitle("Barcode 3 distribution across plate")
    p <- p1 + p2 + p3 + plot_layout(ncol = 2)
    print(p)
  }
  invisible(dev.off())
} 


plot_rna_alignment_stats <- function(input_path, output_path){
  if (!dir.exists(paste0(input_path,"/RNA"))){
    return()
  }
  else{
    rna.align.stats.list <- list.files(paste0(input_path,"/RNA"), pattern="alignment_stats.json", recursive=T)
    rna.align.stats.data <- lapply(paste0(input_path,"/RNA/",rna.align.stats.list),function(f){fromJSON(file=f) %>% as.data.frame()}) %>% do.call(rbind, .) %>% t()
    
    rna.name.list <- unlist(strsplit(paste0("RNA/",rna.align.stats.list), "/alignment_stats.json"))
    rna.name.list <- gsub("/", "_", rna.name.list)
    colnames(rna.align.stats.data) <- rna.name.list
    
    # take the difference
    rows <- c("unmatched_barcode", "unaligned", "multimapped", "unannotated", "duplicated", "post_dedup")
    rna.align.stats.toplot <- data.frame(c(rna.align.stats.data["total_reads",] - rna.align.stats.data["valid_barcode_reads",],
                                           rna.align.stats.data["valid_barcode_reads",] - rna.align.stats.data["total_aligned_reads",],
                                           rna.align.stats.data["multimapped_reads",],
                                           rna.align.stats.data["total_aligned_reads",] - rna.align.stats.data["multimapped_reads",] - rna.align.stats.data["total_annotated_reads",],
                                           rna.align.stats.data["total_annotated_reads",] - rna.align.stats.data["total_deduped_reads",],
                                           rna.align.stats.data["total_deduped_reads",]
    ) )
    colnames(rna.align.stats.toplot) <- "reads"
    
    rna.align.stats.toplot["read_status"] <- factor(lapply(rows,function(n){rep(n,length(rna.name.list))}) %>% unlist, 
                                                    levels = rows)
    rna.align.stats.toplot["sublibrary"] <- rep(rna.name.list, length(rows))
    
    # plot
    p1 <- ggplot(rna.align.stats.toplot, aes(fill=read_status, y=reads, x=sublibrary)) + geom_bar(position="stack",stat="identity") + 
      theme_classic() + theme(text=element_text(family="sans", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      scale_fill_brewer(palette="Blues") 
    p2 <- ggplot(rna.align.stats.toplot, aes(fill=read_status, y=reads, x=sublibrary)) + geom_bar(position="fill",stat="identity") + 
      theme_classic() + theme(text=element_text(family="sans", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      scale_fill_brewer(palette="Blues") + ylab("percentage of reads")
    p <- p1+p2
    
    ggsave(plot=p, paste0(output_path,"/rna_alignment_stats.pdf"), dpi=600, width=10, height=5)
    write.table(rna.align.stats.data, paste0(output_path,"/rna_alignment_stats.tsv"), quote=F, sep="\t")
    return(rna.align.stats.data)
  }
}


plot_atac_alignment_stats <- function(input_path, output_path){
  if (!dir.exists(paste0(input_path,"/ATAC"))){
    return()
  }
  else{
    atac.align.stats.list <- list.files(paste0(input_path,"/ATAC"), pattern="alignment_stats.json", recursive=T)
    atac.align.stats.data <- lapply(paste0(input_path,"/ATAC/",atac.align.stats.list),function(f){fromJSON(file=f) %>% as.data.frame()}) %>% do.call(rbind, .) %>% t()
    
    atac.name.list <- unlist(strsplit(paste0("ATAC/",atac.align.stats.list), "/alignment_stats.json"))
    atac.name.list <- gsub("/", "_", atac.name.list)
    colnames(atac.align.stats.data) <- atac.name.list
    
    # take the difference
    rows <- c("unmatched_barcode", "unaligned", "filtered", "duplicated", "mitochondrial", "postdedup_nomito")
    atac.align.stats.toplot <- data.frame(c(atac.align.stats.data["total_reads",] - atac.align.stats.data["valid_barcode_reads",],
                                            atac.align.stats.data["valid_barcode_reads",] - atac.align.stats.data["total_aligned_reads",],
                                            atac.align.stats.data["total_aligned_reads",] - atac.align.stats.data["total_fragments",],
                                            atac.align.stats.data["total_fragments",] - atac.align.stats.data["unique_fragments",],
                                            atac.align.stats.data["mitochondrial_fragments",],
                                            atac.align.stats.data["unique_fragments",] - atac.align.stats.data["mitochondrial_fragments",]
    ) )
    colnames(atac.align.stats.toplot) <- "reads"
    
    atac.align.stats.toplot["read_status"] <- factor(lapply(rows,function(n){rep(n,length(atac.name.list))}) %>% unlist, 
                                                     levels = rows)
    atac.align.stats.toplot["sublibrary"] <- rep(atac.name.list, length(rows))
  
    # plot
    p1 <- ggplot(atac.align.stats.toplot, aes(fill=read_status, y=reads, x=sublibrary)) + geom_bar(position="stack",stat="identity") + 
          theme_classic() + theme(text=element_text(family="sans", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
          scale_fill_brewer(palette="Blues") 
    p2 <- ggplot(atac.align.stats.toplot, aes(fill=read_status, y=reads, x=sublibrary)) + geom_bar(position="fill",stat="identity") + 
          theme_classic() + theme(text=element_text(family="sans", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
          scale_fill_brewer(palette="Blues") + ylab("percentage of reads")
    p <- p1 + p2
    ggsave(plot=p, paste0(output_path,"/atac_alignment_stats.pdf"), dpi=600, width=10, height=5)
    write.table(atac.align.stats.data, paste0(output_path,"/atac_alignment_stats.tsv"), quote=F, sep="\t")
    return(atac.align.stats.data)
  }
}


plot_total_reads <- function(rna.align.stats.data, atac.align.stats.data, output_path){
  reads <- c(rna.align.stats.data["total_reads",],atac.align.stats.data["total_reads",]) %>% as.data.frame
  colnames(reads) <- "total_reads"
  reads$sublib <- rownames(reads)
  reads <- reads[!(reads$sublib %in% c("ATAC", "RNA")),] # exclude the summary counts of all ATAC/all RNA
  ggplot(reads, aes(y=total_reads, x=sublib)) + geom_bar(position="stack",stat="identity", fill="#094891") + theme_classic() + 
        theme(text=element_text(family="sans", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        ggtitle(paste0("total reads: ", format(sum(reads$total_reads),big.mark=",")))
  
  ggsave(paste0(output_path,"/total_reads.pdf"), width=10, height=5)
  write.table(reads, paste0(output_path, "/total_reads.tsv"), quote=F, sep="\t")
}


########## main ############
bc.stats.list <- list.files(input_path, pattern="barcode_stats.json", recursive=T)
bc.stats.data <-sapply(paste0(input_path,"/",bc.stats.list),function(f){fromJSON(file=f)})
name.list <- unlist(strsplit(bc.stats.list, "/barcode_stats.json"))
name.list <- gsub("/", "_", name.list)
colnames(bc.stats.data) <- name.list

if (!dir.exists(output_path)){
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
}
plot_barcode_stats(bc.stats.data, output_path)
plot_barcode_distr(bc.stats.data, output_path)
rna.align.stats.data <- plot_rna_alignment_stats(input_path, output_path)
atac.align.stats.data <- plot_atac_alignment_stats(input_path, output_path)
plot_total_reads(rna.align.stats.data, atac.align.stats.data, output_path)
