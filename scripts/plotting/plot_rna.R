suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(ggrastr)
  library(gridExtra)
  library(BPCells)
})


########## arguments ############
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)

input_path <- args[1] # path to the shareseq.smk output folder
output_path <- args[2] # path to save script outputs (plots, tables, objects)
flag_saveobj <- args[3] # T or F, whether to save the raw RDS Seurat project

stopifnot(dir.exists(paste0(input_path, "/RNA")))
stopifnot(flag_saveobj %in% c("T","F"))

########## functions ############
read_data <- function(input.dir, use.sublib=FALSE, min.umi=100){
  # read the matrix market files output from https://github.com/GreenleafLab/shareseq-pipeline 
  # input.dir: output folder from the preprocessing pipeline 
  # use.sublib: if true, read from RNA/sublibraries instead of RNA/samples
  # min.umi: filter out cell barcodes with less than min.umi reads
 
  if (use.sublib){
    sample.list <- list.files(file.path(input.dir,"/RNA/sublibraries/"), pattern="matrix.mtx.gz", recursive=T) %>% 
      strsplit("/matrix.mtx.gz") %>% unlist
    
    prefix <- file.path(input.dir,"/RNA/sublibraries/")
    mtx.suffix <- "/matrix.mtx.gz"
    cells.suffix <- "/barcodes.tsv.gz"
    features.suffix <- "/features.tsv.gz"
  } else{
    sample.list <- list.files(file.path(input.dir,"/RNA/samples/"), pattern=".matrix.mtx.gz", recursive=F) %>% 
      strsplit(".matrix.mtx.gz") %>% unlist
    
    prefix <- file.path(input.dir,"/RNA/samples/")
    mtx.suffix <- ".matrix.mtx.gz"
    cells.suffix <- ".barcodes.tsv.gz"
    features.suffix <- ".features.tsv.gz"
  }
  
  data.list <- lapply(sample.list, function(n){
    ReadMtx(mtx=paste0(prefix, n,  mtx.suffix),
            cells=paste0(prefix, n, cells.suffix),
            features=paste0(prefix, n, features.suffix))
  })
  
  names(data.list) <- sample.list
  
  obj.list <- lapply(sample.list,function(n){
    CreateSeuratObject(counts=data.list[[n]], project=n, min.cells=0, min.features=min.umi)
  })
  
  
  if (length(obj.list) == 1){
    proj <- RenameCells(obj.list[[1]], add.cell.id = sample.list)
  }else{
    proj <- merge(obj.list[[1]], y = obj.list[2:length(obj.list)], add.cell.ids = sample.list, project = "proj") 
  }
  
  
  # adding a few attributes
  proj$sublib <- proj$orig.ident
  if (use.sublib){
    tmp <- strsplit(rownames(proj@meta.data), "_")
    proj$sample <- unlist(lapply(tmp, function(n){n[1]}))
    proj$cb <- unlist(lapply(tmp, function(n){n[2]}))
  } else{
    tmp <- strsplit(rownames(proj@meta.data), paste0("_",proj$orig.ident,"_"))
    proj$sample <- unlist(lapply(tmp, function(n){n[1]}))
    proj$cb <- unlist(lapply(tmp, function(n){n[2]}))
  }
  mito_gene_id <- rownames(proj)[grepl("^MT-", rownames(proj))]
  proj$percent.mt <- PercentageFeatureSet(proj, features=mito_gene_id)
  return(proj)
}


plot_knee_all <- function(proj, output_path, umi.cutoff=1000, gene.cutoff=500){
  # overview stats for all cells from all samples
  g1 <- plot_read_count_knee(proj$nCount_RNA, cutoff = umi.cutoff) + 
    ggtitle("nUMI vs cell barcodes, all") + ylab("nUMI") + geom_point(color="darkblue")
  g2 <- plot_read_count_knee(proj$nFeature_RNA, cutoff = gene.cutoff) + 
    ggtitle("nGene vs cell barcodes, all") + ylab("nGene") + geom_point(color="darkblue")
  
  downsample.idx <- plot_read_count_knee(proj$nCount_RNA, return_data=TRUE)
  g3_base <- ggplot(proj@meta.data[names(downsample.idx$data$reads),]) + aes(x=nCount_RNA,y=nFeature_RNA) + 
    geom_point(colour="darkblue") + theme_classic() + 
    xlab("nUMI") + ylab("nGene") + ggtitle("nGene vs nUMI, all") + 
    geom_vline(xintercept=umi.cutoff, linetype="dashed") + geom_hline(yintercept=gene.cutoff, linetype="dashed") +
    scale_x_log10() + scale_y_log10() + annotation_logticks(sides="bl") 
  
  cell_count_text_a <- paste(paste0("Number of Cells that detected > 0 UMIs:       ",sum(proj$nCount_RNA > 0)), 
                             paste0("Number of Cells that detected > 10 UMIs:     ",sum(proj$nCount_RNA > 10)),
                             paste0("Number of Cells that detected > 100 UMIs:   ",sum(proj$nCount_RNA > 100)),
                             paste0("Number of Cells that detected > 500 UMIs:   ",sum(proj$nCount_RNA > 500)),
                             paste0("Number of Cells that detected > 1000 UMIs: ",sum(proj$nCount_RNA > 1000)), sep="\n")
  
  cell_count_text_b <- paste(paste0("Number of Cells that detected > 0 genes:       ",sum(proj$nFeature_RNA > 0)), 
                             paste0("Number of Cells that detected > 10 genes:     ",sum(proj$nFeature_RNA > 10)),
                             paste0("Number of Cells that detected > 100 genes:   ",sum(proj$nFeature_RNA > 100)),
                             paste0("Number of Cells that detected > 500 genes:   ",sum(proj$nFeature_RNA > 500)),
                             paste0("Number of Cells that detected > 1000 genes: ",sum(proj$nFeature_RNA > 1000)), sep="\n")
  
  g3a <- g3_base + labs(caption=cell_count_text_a) + theme(plot.caption = element_text(hjust = 0))
  g3b <- g3_base + labs(caption=cell_count_text_b) + theme(plot.caption = element_text(hjust = 0))
  pdf(paste0(output_path,"/rna_kneeplots_all.pdf"),width=10,height=5)
  print(g1+g2)
  print(g3a+g3b)
  invisible(dev.off())
}


plot_knee_sublib <- function(proj, output_path, umi.cutoff=1000, gene.cutoff=500){
  # UMI and gene knee plots for each sublibrary
  sublib.list <- unique(proj$sublib)
  
  pdf(paste0(output_path, "/rna_kneeplots_persublib.pdf"),width=10,height=5)
  for (sample in sublib.list){
    subproj <- proj@meta.data[proj$sublib == sample,]
    g1 <- plot_read_count_knee(subproj$nCount_RNA, cutoff = umi.cutoff) + 
      ggtitle(paste0("nUMI vs cell barcodes, ", sample)) + ylab("nUMI") + geom_point(color="darkblue")
    g2 <- plot_read_count_knee(subproj$nFeature_RNA, cutoff = gene.cutoff) + 
      ggtitle(paste0("nGene vs cell barcodes, ", sample)) + ylab("nGene") + geom_point(color="darkblue")
    print(g1+g2)
  }
  invisible(dev.off())
}


plot_knee_sample <- function(proj, output_path, umi.cutoff=1000, gene.cutoff=500){
  # UMI and gene knee plots for each sample
  sample.list <- unique(proj$sample)
  
  pdf(paste0(output_path, "/rna_kneeplots_persample.pdf"),width=10,height=5)
  for (sample in sample.list){
    subproj <- proj@meta.data[proj$sample == sample,]
    g1 <- plot_read_count_knee(subproj$nCount_RNA, cutoff = umi.cutoff) + 
      ggtitle(paste0("nUMI vs cell barcodes, ", sample)) + ylab("nUMI") + geom_point(color="darkblue")
    g2 <- plot_read_count_knee(subproj$nFeature_RNA, cutoff = gene.cutoff) + 
      ggtitle(paste0("nGene vs cell barcodes, ", sample)) + ylab("nGene") + geom_point(color="darkblue")
    print(g1+g2)
  }
  invisible(dev.off())
}


plot_violin_umi_gene <- function(proj, output_path){
  median_tbl <- proj@meta.data %>% group_by(sublib) %>% 
        summarise(median_umi=median(nCount_RNA), median_gene=median(nFeature_RNA), cells_passfilter=n())
  p1 <- VlnPlot(proj, features = c("nCount_RNA","nFeature_RNA","percent.mt"), pt.size=0, group.by="sublib", log=T)
  p <- p1 / tableGrob(median_tbl, rows = NULL)
  ggsave(plot=p, file=paste0(output_path,"/rna_umi_gene_violin.pdf"), dpi=300, width=10, height=8)
}


plot_umap <- function(proj, output_path){
  proj <- NormalizeData(proj) %>% FindVariableFeatures(nfeatures=2000) %>% ScaleData %>% RunPCA %>% FindNeighbors(dims = 1:50) %>%
    FindClusters(resolution = 0.2) %>% RunUMAP(dims = 1:50)
  
  p <- DimPlot(proj, group.by="seurat_clusters") +  DimPlot(proj, group.by="sample") +
    FeaturePlot(proj, features="nCount_RNA", min.cutoff="q10", max.cutoff="q90", raster=T, raster.dpi=c(1024,1024)) +
    FeaturePlot(proj, features="TOP2A", min.cutoff="q10", max.cutoff="q90", raster=T, raster.dpi=c(1024,1024))
  ggsave(plot=p, paste0(output_path,"/rna_umap.pdf"), width=10, height=8)
  return(proj)
}

plot_cluster_composition <- function(proj, output_path){
  #colpal <- c("F8766D", "7CAE00","0CB702","00A9FF","FF61CC")
  # per sublib
  g1 <- ggplot(proj@meta.data, aes(x=orig.ident, y=1, fill=seurat_clusters)) + geom_bar(position="fill", stat= "identity") + 
    theme_classic() + theme(text=element_text(family="sans", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    xlab("Sublibrary") + ylab("percentage of cells")
  g1b <- ggplot(proj@meta.data, aes(x=orig.ident, y=1, fill=seurat_clusters)) + geom_bar(position="stack", stat= "identity") + 
    theme_classic() + theme(text=element_text(family="sans", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    xlab("Sublibrary") + ylab("# cells")
  g2 <- ggplot(proj@meta.data, aes(x=seurat_clusters, y=1, fill=orig.ident)) + geom_bar(position="fill", stat= "identity") + 
    theme_classic() + theme(text=element_text(family="sans", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    xlab("cluster") + ylab("percentage of cells")
  g2b <- ggplot(proj@meta.data, aes(x=seurat_clusters, y=1, fill=orig.ident)) + geom_bar(position="stack", stat= "identity") + 
    theme_classic() + theme(text=element_text(family="sans", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    xlab("cluster") + ylab("# cells")
  
  # per sample
  g3 <- ggplot(proj@meta.data, aes(x=sample, y=1, fill=seurat_clusters)) + geom_bar(position="fill", stat= "identity") + 
    theme_classic() + theme(text=element_text(family="sans", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    xlab("Sample") + ylab("percentage of cells")
  g3b <- ggplot(proj@meta.data, aes(x=sample, y=1, fill=seurat_clusters)) + geom_bar(position="stack", stat= "identity") + 
    theme_classic() + theme(text=element_text(family="sans", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    xlab("Sample") + ylab("# cells")
  g4 <- ggplot(proj@meta.data, aes(x=seurat_clusters, y=1, fill=sample)) + geom_bar(position="fill", stat= "identity") + 
    theme_classic() + theme(text=element_text(family="sans", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    xlab("cluster") + ylab("percentage of cells")
  g4b <- ggplot(proj@meta.data, aes(x=seurat_clusters, y=1, fill=sample)) + geom_bar(position="stack", stat= "identity") + 
    theme_classic() + theme(text=element_text(family="sans", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    xlab("cluster") + ylab("# cells")
  
  pdf(paste0(output_path,"/rna_cluster_composition.pdf"), height=5, width=10)
  print(g1+g1b)
  print(g2+g2b)
  print(g3+g3b)
  print(g4+g4b)
  invisible(dev.off())

}


########## main ############
message("reading data")
proj <- read_data(input_path)
#proj <- read_data(input_path, use.sublib=T) # use this if only want sublibrary level qc, no sample info

if (flag_saveobj == "T"){
  message("saving raw object")
  saveRDS(proj, file=paste0(output_path, "/RNA_proj_raw.rds"))
}

if (!dir.exists(output_path)){
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
}

message("plotting knee plots")
plot_knee_all(proj, output_path)
plot_knee_sublib(proj, output_path)
plot_knee_sample(proj, output_path)

# simple filter
proj <- proj[,(proj$nCount_RNA>1000) & (proj$nFeature_RNA>500) & (proj$percent.mt<30)]

plot_violin_umi_gene(proj, output_path)

message("plotting umap")
proj <- plot_umap(proj, output_path)
plot_cluster_composition(proj, output_path)

if (flag_saveobj == "T"){
  message("saving clustered object")
  saveRDS(proj, file=paste0(output_path, "/RNA_proj_clustered.rds"))
}
