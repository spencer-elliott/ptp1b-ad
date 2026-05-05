#
# Analyzing data with Seurat
#

# Seurat package contains tools for analyzing scRNA-seq data
library(Seurat)
# Matrix package contains tools for working with matrix classes
library(Matrix)

# Transposing the raw matrix
counts <- t(readMM("counts.mtx"))

# Reading in the barcodes and genes
barcodes <- readLines("counts.barcodes.txt")
features <- readLines("counts.genes.txt")

# Removing version numbers from gene IDs
features <- sub("\\.\\d+$", "", features)

# Assigning row and column names
rownames(counts) <- features
colnames(counts) <- barcodes

# Filtering the barcodes with at least 500 UMIs
counts_filtered <- counts[, colSums(counts) >= 500]

# Creating a Seurat object
so <- CreateSeuratObject(counts = counts_filtered,
                         project = "my_project",
                         min.cells = 3,
                         min.features = 200)
