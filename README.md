# rnaseq_pipeline

This is a complete pipeline for analysing scRNA-seq datasets from 10x Genomics.

The workflow uses Seurat (v3.0.2) and performs all necessary steps for an basic exploratory analysis of your data. These include:
- normalization and scaling of gene-expression matrix, 
- filtering of low quality cells and outliers, 
- dimensionality reduction (PCA and TSNE), 
- clustering (using different resolutions), 
- and differential expression analysis.

Seurat offers a wide variety of functions and tutorials that allow you to perform exhaustive exploration of your scRNA-seq dataset. For this reason, the present pipeline is essentially following the same workflow provided by Satija lab ([Seurat - Getting Started](https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html))

The advantage of this pipeline is that it allows for automatised exploratory analysis of your data, using a one-liner on the CLI. 

```
Rscript script.R -i path/to/input -o path/to/output 
```

In order to use this pipeline you will need to download the repo to your local and set your working directory.
- Click on **Clone or Download**
- Download ZIP
- Change working directory

```
cd /path/to/repo 
```

## Prerequisites
- Install the latest [R](https://www.r-project.org/) (I have used R v3.6.0)
- Packages needed for running the pipeline will automatically be installed, if not already present. If Seurat is already installed, make sure to upgrate to version >= 3.0.0.
- Input option requires the gene-expression files from CellRanger to be unzipped (e.g., `gunzip path/to/files/*.gz` and named as: **barcodes.tsv, genes.tsv, matrix.mtx**

