#!/usr/bin/env Rscript

# Load packages & color palettes
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(viridis))

my.viridis = c(alpha("grey85", alpha = 0.4),
               viridis(10, alpha = 1, begin = 0))

# Parse arguments
option_list = list(
  make_option(c("-i", "--input"), default = "path"),
  make_option(c("-o", "--output"), default = "directory")
  )

options = parse_args(OptionParser(option_list=option_list))

gbm=Read10X(options$input) # Load gbm

mito.genes <- grep(pattern = "^MT-", x = rownames(gbm), 
                   value = F, ignore.case = T)
percent.mito <- Matrix::colSums(gbm[mito.genes,])/Matrix::colSums(gbm)

# Create Seurat object, and add percent.mito to object@meta.data in the percent.mito column
obj <- CreateSeuratObject(raw.data = gbm, 
                          meta.data = data.frame(percent.mito = percent.mito))

dir.create(options$output)
setwd(options$output)

pdf("prefilter_vln_nGene_nUMI_mito.pdf", width = 6, height = 4)
VlnPlot(obj, c("nGene", "nUMI", "percent.mito"),
        nCol = 3, group.by = "orig.ident",
        point.size.use = 0.02, size.x.use = 0)
dev.off()

saveRDS(obj, "prefilter.rds")

print(paste("Filtered", 
            length(which(obj@meta.data$percent.mito > 0.2)), 
            "cells with mitochondrial content > 20%"))

# filter cells with percent.mito >20%
obj <- FilterCells(obj, subset.names = "percent.mito", 
                   low.thresholds = -Inf, high.thresholds = 0.2)

# filter possible doublets based on gene number
outlier.ngene <- mean(obj@meta.data$nGene) + 3*sd(obj@meta.data$nGene)
obj <- FilterCells(obj, subset.names = "nGene", 
                   low.thresholds = -Inf, high.thresholds = outlier.ngene)
print(paste("Filtered",
            length(which(obj@meta.data$nGene > outlier.ngene)),
            "cells with nGene > 3*sd"))

# plot after filtering of low quality cells and doublets
pdf("postfilter_vln_nGene_nUMI_mito.pdf", width = 6, height = 4)
VlnPlot(obj, c("nGene", "nUMI", "percent.mito"),
        nCol = 3, group.by = "orig.ident",
        point.size.use = 0.02, size.x.use = 0)
dev.off()

saveRDS(obj, "filtered.rds")

# Normalize and scale
obj <- NormalizeData(object = obj,
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4)
obj <- ScaleData(object = obj, vars.to.regress = c("nUMI", "percent.mito"))

obj <- FindVariableGenes(obj, selection.method = "dispersion", top.genes = 2000)

saveRDS(obj, "scaled.rds")

# Run PCA and tSNE for clustering cell types
obj <- RunPCA(obj, pc.genes = obj@var.genes, pcs.compute = 30)

pc_sdev <- GetDimReduction(obj, slot = "sdev")

obj <- RunTSNE(obj, dims.use = 1:max(which(pc_sdev >= 2)))

pdf("tsne.pdf", width = 6, height = 6)
TSNEPlot(obj, no.legend = T,
         plot.title = paste("dim.use: top", max(which(pc_sdev >= 2)), sep = " "))
dev.off()

saveRDS(obj, "tsne.rds")

obj <- RunUMAP(obj, dims.use = 1:max(which(pc_sdev >= 2)))
DimPlot(obj, no.legend = T, reduction.use = "umap",
         plot.title = paste("dim.use: top", max(which(pc_sdev >= 2)), sep = " "))
