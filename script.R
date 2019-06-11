#!/usr/bin/env Rscript

# Load packages & color palettes
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(viridis))

my.viridis = c(alpha("grey85", alpha = 0.4),
               viridis(10, alpha = 1, begin = 0))

option_list = list(
  make_option(c("-i", "--input"), default="path")
  )

options = parse_args(OptionParser(option_list=option_list))
gbm=Read10X(options$input)

print(gbm[1:3,1:3])

# Load gbm

# 