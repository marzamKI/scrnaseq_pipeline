#!/usr/bin/env Rscript
options(warn=-1)

if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages("Seurat")
if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse")
if (!requireNamespace("cowplot", quietly = TRUE))
  install.packages("cowplot")

suppressPackageStartupMessages(require(optparse))

# Parse arguments
option_list = list(
  make_option(c("-i", "--input"), default = "path",
              help = "Set the path to the directory containing
              matrix.mtx, genes.tsv, and barcodes.tsv files. 
              Make sure files are unzipped"),
  make_option(c("-o", "--output"), default = "directory",
              help = "Specify output directory for storing plots and tables"),
  make_option(c("-g", "--genes"), default = NA, type = "character",
              help = "(Optional) Specify path to a .csv file listing genes"),
  make_option(c("-c", "--cellcycle"), default = NA, type = "character",
              help = "(Optional) If TRUE, calculates cell cycle score using 
              S- and G2M-related genes")
)

options = parse_args(OptionParser(option_list=option_list), positional_arguments = F)


# Load packages & color palettes
print("Loading packages")
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(cowplot))

print("Loading gene expression matrix and creating Seurat object")
gbm=Read10X(options$input) # Load gbm

# Create Seurat object, and add percent.mito to object@meta.data in the percent.mito column
obj <- CreateSeuratObject(gbm, min.features = 200)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

dir.create(options$output)
setwd(options$output)
print(paste("Saving output in ", getwd(), sep = ""))

prefilter <- VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     ncol = 3, pt.size = 0.2) 

saveRDS(obj, "prefilter.rds")

print("Filter cells with percent.mito >20% and nfeatures > 3sd from the mean")
print(paste("Filtered", 
            length(which(obj@meta.data$percent.mt > 20)), 
            "cells with mitochondrial content > 20%")
)

outlier.ngene <- mean(obj@meta.data$nFeature_RNA) + 3*sd(obj@meta.data$nFeature_RNA)
print(paste("Filtered",
            length(which(obj@meta.data$nFeature_RNA > outlier.ngene)),
            "cells with nGene > 3*sd")
)

obj <- subset(obj, subset = nFeature_RNA < outlier.ngene & percent.mt < 20)

postfilter <- VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                      ncol = 3, pt.size = 0.2) 

# plot after filtering of low quality cells and doublets
pdf("PrePostFilterVln.pdf", width = 6, height = 9)
plot_grid(prefilter, postfilter, ncol = 1)
dev.off()

saveRDS(obj, "filtered.rds")

print(paste("Total number of cells to be analysed:",
      ncol(obj@assays$RNA)
      )
)

# Normalize and scale
print("Log-normalizing data and looking for HVGs")
obj <- NormalizeData(object = obj,
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4)

obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, ynudge = 0, xnudge = 0)

pdf("HVG.pdf", width = 6, height = 4)
plot2
dev.off()

print("Scaling data using SCTransform")
obj <- SCTransform(obj, vars.to.regress = c("nFeature_RNA","percent.mt"), verbose = FALSE)

saveRDS(obj, "scaled.rds")

# Run PCA
print("Performing dimensionality reduction")
obj <- RunPCA(obj, features = VariableFeatures(object = obj))

write.csv(obj@reductions$pca@feature.loadings, "PCAFeatureLoadings.csv")
write.csv(obj@reductions$pca@cell.embeddings, "PCACellEmbeddings.csv")

obj <- JackStraw(obj, num.replicate = 100)
obj <- ScoreJackStraw(obj, dims = 1:20)

print("Jackstraw was used to calculate p-value of the top 20 PCs. See below")
JS_pvalue <- as.data.frame(JS(object = obj[['pca']], slot = 'overall'))
print(JS_pvalue)

# Run clustering analysis
print("Significant PCs are used for subsequent tSNE and clustering analyses")
sig_dims <- which(JS_pvalue$Score < 0.05)
obj <- FindNeighbors(obj, dims = sig_dims)

print("Running clustering analysis at resolutions 0.1 to 1")
# Run clustering using different resolutions
for(res in c((1:10)/10)){
  obj <- FindClusters(obj, resolution = res)
}

# Run tSNE
obj <- RunTSNE(obj, dims = sig_dims, check_duplicates = F)

saveRDS(obj, "tsne.rds")

resolutions <- paste("SCT_snn_res.", (1:10)/10, sep="")
res_plots <- list()
for(res in resolutions){
  res_plots[[res]]<-DimPlot(obj, reduction = "tsne", label = T, label.size = 6, group.by = res) +
    ggtitle(paste("dim.use:", max(which(JS_pvalue$Score < 0.05)), 
                  "res:", res, sep = " ")) +
    theme(legend.position = "none")
}

pdf("tSNEPlot.pdf", width = 12, height = 16)
plot_grid(plotlist = res_plots, ncol = 3)
dev.off()

write.csv(obj@reductions$tsne@cell.embeddings, "tSNECellEmbeddings.csv")

# Run UMAP
obj <- RunUMAP(obj, dims = sig_dims)

saveRDS(obj, "umap.rds")

resolutions <- paste("SCT_snn_res.", (1:10)/10, sep="")
res_plots <- list()
for(res in resolutions){
  res_plots[[res]]<-DimPlot(obj, reduction = "umap", label = T, label.size = 6, group.by = res) +
    ggtitle(paste("dim.use:", max(which(JS_pvalue$Score < 0.05)), 
                  "res:", res, sep = " ")) +
    theme(legend.position = "none")
}

pdf("UMAPPlot.pdf", width = 12, height = 16)
plot_grid(plotlist = res_plots, ncol = 3)
dev.off()

write.csv(obj@reductions$umap@cell.embeddings, "UMAPCellEmbeddings.csv")


if (!is.na(options$genes)) {
  file <- read.csv("genes.csv", header = T, stringsAsFactors = F)
  genes <- file[,1]
  pdf("FeaturePlot.pdf", paper = 'special')
  for (i in 1:length(x = genes)) {
    print(FeaturePlot(obj, reduction = "umap", features = genes[[i]]))}
  dev.off()
} else {
  print("No list of genes was provided for FeaturePlot() function")
}

### Cell cycle

if(!is.na(options$cellcycle)) { 
  print("Scoring cell cycle module activity")
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  pdf("CellCycleScore.pdf", paper = 'special', width = 9, height = 4)
  print(FeaturePlot(obj, reduction = "umap", features = c("S.Score", "G2M.Score")))
  dev.off()
}

