
# ==========
# PTP1B scRNA-seq Analysis
# (One sample)
# By Spencer :D
# ==========


# Setup ----

# Load packages
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(here)

# Set the working directory to Kallisto output folder
setwd(here("results/kallisto_output"))

# Save random seed for rerunning
seed <- readRDS(here("results/original_seed.rds"))
assign(".Random.seed", seed, envir = .GlobalEnv)


# Creating the Seurat object ----

# Read in the raw matrix
counts <- readMM("counts.mtx")

# Transpose the matrix
counts <- t(counts)

# Read in the barcodes and genes
barcodes <- readLines("counts.barcodes.txt")
features <- readLines("counts.genes.txt")

# Remove version numbers from gene IDs
features <- sub("\\.\\d+$", "", features)

# Assign row and column names
rownames(counts) <- features
colnames(counts) <- barcodes

# Filter the barcodes
counts_filtered <- counts[, colSums(counts) >= 500]

# Create the Seurat object
so <- CreateSeuratObject(counts = counts_filtered,
                         project = "ptp1b_ad",
                         min.cells = 3,
                         min.features = 200)


# QC / Filtering ----

# Read in the annotation file
gtf <- read.table(here("results/gencode.vM23.annotation.gtf"),
                  sep = "\t",
                  comment.char = "#",
                  header = FALSE)

# Create a list of mitochondrial gene IDs
mt.gtf <- gtf[gtf$V1 == "chrM" & gtf$V3 == "gene", ]
attr <- mt.gtf$V9
mt.genes <- sub('.*gene_id "?([^";]+)"?.*', '\\1', attr)
mt.genes <- sub("\\..*", "", mt.genes)

# Calculate mitochondrial percentages
so[["percent.mt"]] <- PercentageFeatureSet(so, features = mt.genes)

# Visualize QC with a violin plot
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualize QC with a feature-scatter plot
plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter the low-quality cells
so <- subset(so, subset = nFeature_RNA > 200 &
               nFeature_RNA < 5000 &
               percent.mt < 15)


# Preprocessing ----

# Normalize the data
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)

# Map gene names to their Ensembl IDs
gene.gtf <- gtf[gtf$V3 == "gene",]
gene.ids   <- sub('.*gene_id "?([^";]+)"?.*',   '\\1', gene.gtf$V9)
gene.ids <- sub("\\.\\d+$", "", gene.ids)
gene.names <- sub('.*gene_name "?([^";]+)"?.*', '\\1', gene.gtf$V9)
gene.map <- data.frame(ensembl_id = gene.ids,
                       gene_name  = gene.names)

# Rename all Ensembl IDs to gene names
new.names <- gene.map$gene_name[match(rownames(so), gene.map$ensembl_id)]
rownames(so[["RNA"]]) <- make.unique(new.names)

# Find the 2000 most variable genes
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
# Find the top 10 most variable genes
top10 <- head(VariableFeatures(so), 10)

# Visualize the results as a variable feature plot
plot1 <- VariableFeaturePlot(so)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Scale the data
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)

# Run the PCA
so <- RunPCA(so, features = VariableFeatures(so))

# Visualize the PCs with an elbow plot
ElbowPlot(so)


# Clustering ----

# Cluster the cells
so <- FindNeighbors(so, dims = 1:15)
so <- FindClusters(so, resolution = 0.5)

# Condense all the PCs down to 2 dimensions (UMAP)
so <- RunUMAP(so, dims = 1:15)

# Visualize the clusters with a dimensional reduction plot
DimPlot(so, reduction = "umap", label = TRUE)

# Visualize expression of one gene with a feature plot
FeaturePlot(so, features = "Gsx1")

# Identify the marker genes
so.markers <- FindAllMarkers(so, only.pos = TRUE)

# List the top 5 marker genes for each cluster
top5 <- so.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Visualize marker genes with a heatmap
DoHeatmap(so, features = top5$gene) + NoLegend()


# Automatic annotations (SingleR) ----

# Libraries for SingleR
library(SingleR)
library(celldex)

# Reference dataset for SingleR
ref <- celldex::MouseRNAseqData()

# Make the vector of cell type predictions
counts.singleR <- GetAssayData(so, layer = "counts")
pred <- SingleR(test = counts.singleR, ref = ref, labels = ref$label.main)
head(pred$labels, 10)

# Add the predictions to the Seurat object, and plot with these labels
so$SingleR.labels <- pred$labels
DimPlot(so, group.by = "SingleR.labels", reduction = "umap", label = TRUE)


# Previewing next steps ----

# Viewing Ptpn1 expression across cell types
FeaturePlot(so, features = "Ptpn1")

# Snap25 (neuron marker gene)
FeaturePlot(so, features = "Snap25")

