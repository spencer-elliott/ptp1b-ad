#
# Analyzing data with Seurat
#

# Setting the working directory to Kallisto output folder
setwd("results/kallisto_output")

# Loading packages
library(Seurat)
library(Matrix)
library(dplyr)
library(SingleR)
library(celldex)

# Reading in the raw matrix
counts <- readMM("counts.mtx")

# Transposing the matrix
counts <- t(counts)

# Reading in the barcodes and genes
barcodes <- readLines("counts.barcodes.txt")
features <- readLines("counts.genes.txt")

# Removing version numbers from gene IDs
features <- sub("\\.\\d+$", "", features)

# Assigning row and column names
rownames(counts) <- features
colnames(counts) <- barcodes

# Filtering the barcodes
counts_filtered <- counts[, colSums(counts) >= 500]

# Creating the Seurat object
so <- CreateSeuratObject(counts = counts_filtered,
                         project = "ptp1b_ad",
                         min.cells = 3,
                         min.features = 200)

# Reading in the annotation file
gtf <- read.table("../gencode.vM23.annotation.gtf",
                  sep = "\t",
                  comment.char = "#",
                  header = FALSE)

# Creating a list of mitochondrial gene IDs
mt.gtf <- gtf[gtf$V1 == "chrM" & gtf$V3 == "gene", ]
attr <- mt.gtf$V9
mt.genes <- sub('.*gene_id "?([^";]+)"?.*', '\\1', attr)
mt.genes <- sub("\\..*", "", mt.genes)

# Calculating the mitochondrial percentage
so[["percent.mt"]] <- PercentageFeatureSet(so, features = mt.genes)

# Visualizing QC with a violin plot
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualizing QC with a feature-scatter plot
plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filtering out low-quality cells
so <- subset(so, subset = nFeature_RNA > 200 &
               nFeature_RNA < 5000 &
               percent.mt < 15)

# Normalizing the data
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)

# Mapping gene names to the Ensembl IDs
gene.gtf <- gtf[gtf$V3 == "gene",]
gene.ids   <- sub('.*gene_id "?([^";]+)"?.*',   '\\1', gene.gtf$V9)
gene.ids <- sub("\\.\\d+$", "", gene.ids)
gene.names <- sub('.*gene_name "?([^";]+)"?.*', '\\1', gene.gtf$V9)
gene.map <- data.frame(ensembl_id = gene.ids,
                       gene_name  = gene.names)

# Renaming all Ensembl IDs to gene names
new.names <- gene.map$gene_name[match(rownames(so), gene.map$ensembl_id)]
rownames(so[["RNA"]]) <- make.unique(new.names)

# Finding the 2000 most variable genes
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
# Finding the top 10 most variable genes
top10 <- head(VariableFeatures(so), 10)

# Visualizing the results as a variable feature plot
plot1 <- VariableFeaturePlot(so)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Scaling the data
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)

# Running the PCA
so <- RunPCA(so, features = VariableFeatures(so))

# Visualizing the PCs with an elbow plot
ElbowPlot(so)

# Clustering the cells
so <- FindNeighbors(so, dims = 1:15)
so <- FindClusters(so, resolution = 0.5)

# Condensing all the PCs down to 2 dimensions
so <- RunUMAP(so, dims = 1:15)

# Visualizing the clusters with a dimensional reduction plot
DimPlot(so, reduction = "umap", label = TRUE)

# Visualizing expression of one gene with a feature plot
FeaturePlot(so, features = "Gsx1")

# Identifying the marker genes
so.markers <- FindAllMarkers(so, only.pos = TRUE)

# Listing the top 5 marker genes for each cluster
top5 <- so.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Visualizing marker genes with a heatmap
DoHeatmap(so, features = top5$gene) + NoLegend()
