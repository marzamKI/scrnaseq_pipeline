#!/usr/bin/env Rscript
options(warn=-1)

suppressPackageStartupMessages(require(optparse))

# Parse arguments
option_list = list(
  make_option(c("-i", "--input"), default = "path",
              help = "Set the path to the R objects previously saved."),
  make_option(c("-o", "--output"), default = NA, type = "character",
              help = "(Optional) Specify output directory for storing plots and tables.")
)

options = parse_args(OptionParser(option_list=option_list), positional_arguments = F)

# Load packages & color palettes
print("Loading packages")
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(cowplot))
suppressPackageStartupMessages(require(tidyverse))

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

