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

obj <- CreateSeuratObject(gbm) # Create object

# mitochondrial genes start with 'mt-'
mito.genes <- grep(pattern = "^MT-", x = rownames(x = obj@data), 
                   value = TRUE, ignore.case = T)
percent.mito <- Matrix::colSums(obj@raw.data[mito.genes,])/Matrix::colSums(obj@raw.data)
obj_prefilter <- AddMetaData(obj, metadata = percent.mito, col.name = "percent.mito")

dir.create(options$output)
setwd(options$output)

pdf("prefilter_vln_nGene_nUMI_mito.pdf", width = 6, height = 4)
VlnPlot(obj_prefilter, c("nGene", "nUMI", "percent.mito"),
        nCol = 3, group.by = "orig.ident",
        point.size.use = 0.02, size.x.use = 0)
dev.off()

# filter cells with percent.mito >20%
obj <- FilterCells(obj_prefilter, subset.names = "percent.mito", 
                   low.thresholds = -Inf, high.thresholds = 0.2)

# filter possible doublets
outlier.ngene <- mean(obj@meta.data$nGene) + 3*sd(obj@meta.data$nGene)
obj <- FilterCells(obj, subset.names = "nGene", 
                   low.thresholds = -Inf, high.thresholds = outlier.ngene)

# plot after filtering of low quality cells and doublets
pdf("postfilter_vln_nGene_nUMI_mito.pdf", width = 6, height = 4)
VlnPlot(obj, c("nGene", "nUMI", "percent.mito"),
        nCol = 3, group.by = "orig.ident",
        point.size.use = 0.02, size.x.use = 0)
dev.off()

