# scrnaseq_pipeline

This is a complete pipeline for analysing scRNA-seq datasets from 10x Genomics.

The workflow uses Seurat (v3.0.2) and performs all necessary steps for an basic exploratory analysis of your data. These include:
- normalization and scaling of gene-expression matrix, 
- filtering of low quality cells and outliers, 
- dimensionality reduction (PCA, TSNE, and UMAP), 
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
## Overview
The pipeline is composed of two executable R scripts:
- **explore.R** Performs all the exploratory analysis (from loading the gene-expression matrix to UMAP and clustering).
- **cluster.R** Performs differential expression analysis based on the selected resolution used for clustering.

Run `Rscript explore.R -h` to see arguments to be passed with the **explore.R** script:

```
Usage: explore.R [options]


Options:
	-i INPUT, --input=INPUT
		Set the path to the directory containing
              matrix.mtx, genes.tsv, and barcodes.tsv files. 
              Make sure files are unzipped

	-o OUTPUT, --output=OUTPUT
		Specify output directory for storing plots and tables

	-g GENES, --genes=GENES
		(Optional) Specify path to a .csv file listing genes

	-c CELLCYCLE, --cellcycle=CELLCYCLE
		(Optional) If TRUE, calculates cell cycle score using 
              S- and G2M-related genes

	-h, --help
		Show this help message and exit

```
Specify `-o path/to/out` to create an `out` directory where all plots will be saved.
When `-g genes.csv` is added, the script returns a .pdf file displaying gene expression levels on a UMAP embedding. Gene lists should be saved in the repo directory, as a one-column .csv file. E.g.:

```
gene
ASCL1
DCX
MKI67
```

Run `Rscript cluster.R -h` to see arguments to be passed with the **cluster.R** script:

```
Usage: cluster.R [options]


Options:
	-i INPUT, --input=INPUT
		Set the path to the R objects previously saved.

	-r RESOLUTION, --resolution=RESOLUTION
		Set the resolution for the clustering algorithm (i.e., from 0.1 to 1

	-o OUTPUT, --output=OUTPUT
		(Optional) Specify output directory for storing plots and tables.

	-g GENES, --genes=GENES
		(Optional) Specify path to a .csv file listing genes 

	-h, --help
		Show this help message and exit

```
Use `-o`to save output files from **cluster.R** in a new directory. If not specified, files will be saved in the `
out` directory created while running **explore.R**.

## Prerequisites
- Install the latest [R](https://www.r-project.org/) (I have used R v3.6.0)
- Packages needed for running the pipeline will automatically be installed, if not already present. If Seurat is already installed, make sure to upgrate to version >= 3.0.0.
- Input option requires the gene-expression files from CellRanger to be unzipped (e.g., `gunzip path/to/files/*.gz` and named as: **barcodes.tsv, genes.tsv, matrix.mtx**.

## Output 
Output files for **explore.R** include:
- Exploratory plots: 
	- **PrePostFilterVln.pdf**
	- **HVG.pdf**
	- **tSNEPlot.pdf**
	- **UMAPPlot.pdf**
	- **CellCycleScore.pdf**
- Exploratory tables:
	- **PCACellEmbeddings.csv**
	- **PCAFeatureLoadings.csv**
	- **tSNECellEmbeddings.csv**
	- **UMAPCellEmbeddings.csv**
- Seurat objects:
	- **prefilter.rds**
	- **filtered.rds**
	- **scaled.rds**
	- **tsne.rds**
- (Optional) Feature plots (**FeaturePlot.pdf**).

Output files for **cluster.R** include:
- (Optional) Feature plots (**FeaturePlot.pdf**)
- List of cluster markers (**ClusterMarkers.csv**)
- Top 20 cluster markers (**Top20ClusterMarkers.csv**)
- Heatmap with the top 20 features per cluster (**Top20MarkerHM.pdf**)

## Progress
Use ` > progress.log 2>&1` to build a progress report. E.g., 

```
Rscript cluster.R -i path/to/input -g genes.csv > progress.log 2>&1

