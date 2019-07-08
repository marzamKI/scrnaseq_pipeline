#!/usr/bin/env Rscript
options(warn=-1)

suppressPackageStartupMessages(require(optparse))

# Parse arguments
option_list = list(
  make_option(c("-i", "--input"), default = "path",
              help = "Set the path to the R objects previously saved."),
  make_option(c("-r", "--resolution"), default = NA, type = "character",
              help = "Set the resolution for the clustering algorithm (i.e., from 0.1 to 1"),
  make_option(c("-o", "--output"), default = NA, type = "character",
              help = "(Optional) Specify output directory for storing plots and tables."),
  make_option(c("-g", "--genes"), default = NA, type = "character",
              help = "(Optional) Specify path to a .csv file listing genes ")
)

options = parse_args(OptionParser(option_list=option_list), positional_arguments = F)

# Load packages & color palettes
print("Loading packages")
if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(cowplot))
suppressPackageStartupMessages(require(tidyverse))


if (!is.na(options$genes)) {
  file <- read.csv("genes.csv", header = T, stringsAsFactors = F)
  genes <- file[,1]
} else {
  print("No list of genes was provided for FeaturePlot() function")
}

input_dir <- options$input
obj=readRDS(paste(input_dir, "tsne.rds", sep = "")) # Load rds

if (!is.na(options$output)) {
  dir.create(options$output)
  setwd(options$output)
} else {
  setwd(options$input)
}

print(paste("Saving data in",
            getwd(), sep = " "))

if (!is.na(options$genes)) {
  pdf("FeaturePlot.pdf", paper = 'special')
  features <- genes[[1]]
  for (i in 1:length(x = genes)) {
    print(FeaturePlot(obj, features = genes[[i]]))}
  dev.off()
}

# Differential expression analysis at specified clustering resolution
if (!is.na(options$resolution)) {
  Idents(obj) = paste("SCT_snn_res.", options$resolution, sep="")
  obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(obj.markers, "ClusterMarkers.csv")
  top20 <- obj.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
  write.csv(top20, "Top20ClusterMarkers.csv")
  pdf("Top20MarkerHM.pdf", paper = 'special', height = 10, width = 10)
  print(DoHeatmap(obj, features = top20$gene) + NoLegend())
  dev.off()
} else {
  print("Provide clustering resolution to run differential expression analysis")
}

