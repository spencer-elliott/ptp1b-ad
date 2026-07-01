#
# Analyzing data with Seurat
# All four samples
#

#
# Annotating the cell clusters
#

# Setting the working directory to CellRanger output folder
setwd("~/ptp1b_ad/results/cellranger_output")

# Loading packages
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

# Read in the raw matrices
KO1 <- Read10X("KO_SRR36691286/filtered_feature_bc_matrix/")
KO2 <- Read10X("KO_SRR36691287/filtered_feature_bc_matrix/")
AD1 <- Read10X("AD_SRR36691288/filtered_feature_bc_matrix/")
AD2 <- Read10X("AD_SRR36691289/filtered_feature_bc_matrix/")

# Random seed
seed <- readRDS("~/ptp1b_ad/results/original_seed.rds")
assign(".Random.seed", seed, envir = .GlobalEnv)

# Creating a Seurat object for each sample
KO1 <- CreateSeuratObject(counts = KO1, project = "KO", min.cells = 3, min.features = 200)
KO2 <- CreateSeuratObject(counts = KO2, project = "KO", min.cells = 3, min.features = 200)
AD1 <- CreateSeuratObject(counts = AD1, project = "AD", min.cells = 3, min.features = 200)
AD2 <- CreateSeuratObject(counts = AD2, project = "AD", min.cells = 3, min.features = 200)

# Merging the Seurat objects into one
so <- merge(KO1, y = list(KO2, AD1, AD2), 
            add.cell.ids = c("KO1", "KO2", "AD1", "AD2"))

# Tracking mitochondrial percentage
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")

# Visualizing QC
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
# Plot # of UMIs vs mitochondrial content
plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt")
# Plot # of UMIs vs # of features
plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filtering the low quality cells
so <- subset(so, subset = nFeature_RNA > 500 & 
                       nFeature_RNA < 6000 & 
                       nCount_RNA > 1000 &
                       nCount_RNA < 30000 & 
                       percent.mt < 20)

# Normalization
so <- NormalizeData(so)

# Finding the 2000 most variable features
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
# Finding the top 10 variable features
top10 <- head(VariableFeatures(so), 10)

# Plotting the variable features with top 10 labeled
plot1 <- VariableFeaturePlot(so)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Scaling the data
so <- ScaleData(so)

# Running the PCA
so <- RunPCA(so)

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
FeaturePlot(so, features = "Ptpn1")

# Identifying the marker genes (must join layers first)
so <- JoinLayers(so)
so.markers <- FindAllMarkers(so, only.pos = TRUE)

# Listing the top 5 marker genes for each cluster
top5 <- so.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Visualizing marker genes with a heatmap
DoHeatmap(so, features = top5$gene) + NoLegend()

# Libraries for SingleR
library(SingleR)
library(celldex)

# Reference dataset for SingleR
ref <- celldex::MouseRNAseqData()
# Try brain-specific as well?
# ref <- celldex::AllenReferenceBrain()

# Making the vector of cell type predictions
counts.singleR <- GetAssayData(so, layer = "counts")
pred <- SingleR(test = counts.singleR, ref = ref, labels = ref$label.main)
head(pred$labels, 10)

# Adding the predictions to the Seurat object, and plotting again
so$SingleR.labels <- pred$labels
DimPlot(so, group.by = "SingleR.labels", reduction = "umap", label = TRUE)

# Plotting features in the different groups
FeaturePlot(so, features = "Ptpn1", split.by = "orig.ident")
FeaturePlot(so, features = "Syk", split.by = "orig.ident")


# Using the same marker genes as the supplemental info to identify cell types
# (Fig.S4B)

# Astrocyte: Aqp4, Acsbg1, Slc1a2
FeaturePlot(so, features = "Aqp4", label = TRUE)
FeaturePlot(so, features = "Acsbg1", label = TRUE)
FeaturePlot(so, features = "Slc1a2", label = TRUE)
# Clusters: 5, 17, 22

# Choroid plexus cell: Ttr, Ecrg4, Enpp2
FeaturePlot(so, features = "Ttr", label = TRUE)
FeaturePlot(so, features = "Ecrg4", label = TRUE)
FeaturePlot(so, features = "Enpp2", label = TRUE)
# Clusters: 0, 18

# Endothelial cell: Flt1, Slco1a4, Mecom
FeaturePlot(so, features = "Flt1", label = TRUE)
FeaturePlot(so, features = "Slco1a4", label = TRUE)
FeaturePlot(so, features = "Mecom", label = TRUE)
# Clusters: 13

# Macrophage: F13a1, H2-Aa, Lyve1
FeaturePlot(so, features = "F13a1", label = TRUE)
FeaturePlot(so, features = "H2-Aa", label = TRUE)
FeaturePlot(so, features = "Lyve1", label = TRUE)
# Clusters: 11, 20

# Microglia: P2ry12, Hexb, Trem2
FeaturePlot(so, features = "P2ry12", label = TRUE)
FeaturePlot(so, features = "Hexb", label = TRUE)
FeaturePlot(so, features = "Trem2", label = TRUE)
# Clusters: 1, 4, 8, 14, 19

# Neuron: Syt1, Snap25, Nrn1
FeaturePlot(so, features = "Syt1", label = TRUE)
FeaturePlot(so, features = "Snap25", label = TRUE)
FeaturePlot(so, features = "Nrn1", label = TRUE)
# Clusters: 10

# Oligodendrocyte: Ermn, Mog, Aspa
FeaturePlot(so, features = "Ermn", label = TRUE)
FeaturePlot(so, features = "Mog", label = TRUE)
FeaturePlot(so, features = "Aspa", label = TRUE)
# Clusters: 2, 3, 6, 12, 15, 21, 26

# OPC: Pdgfra, Cacng4, Vcan
FeaturePlot(so, features = "Pdgfra", label = TRUE)
FeaturePlot(so, features = "Cacng4", label = TRUE)
FeaturePlot(so, features = "Vcan", label = TRUE)
# Clusters: 9

# T cell: Cd3d, Cd3e, Cd3g
FeaturePlot(so, features = "Cd3d", label = TRUE)
FeaturePlot(so, features = "Cd3e", label = TRUE)
FeaturePlot(so, features = "Cd3g", label = TRUE)
# Clusters: 7, 23, 24, 25

# Cluster 16 - Nothing manual - SingleR says fibroblasts? Not mentioned in paper, drop?

# Cluster annotations based on these results
cluster_annotations <- c(
  "0" = "Choroid plexus cell",
  "1" = "Microglia",
  "2" = "Oligodendrocyte",
  "3" = "Oligodendrocyte",
  "4" = "Microglia",
  "5" = "Astrocyte",
  "6" = "Oligodendrocyte",
  "7" = "T cell",
  "8" = "Microglia",
  "9" = "OPC",
  "10" = "Neuron",
  "11" = "Macrophage",
  "12" = "Oligodendrocyte",
  "13" = "Endothelial cell",
  "14" = "Microglia",
  "15" = "Oligodendrocyte",
  "16" = "Fibroblast",
  "17" = "Astrocyte",
  "18" = "Choroid plexus cell",
  "19" = "Microglia",
  "20" = "Macrophage",
  "21" = "Oligodendrocyte",
  "22" = "Astrocyte",
  "23" = "T cell",
  "24" = "T cell",
  "25" = "T cell",
  "26" = "Oligodendrocyte"
)
so@meta.data$cell_type <- cluster_annotations[as.character(so$seurat_clusters)]
Idents(so) <- "cell_type"
DimPlot(so, reduction = "umap", label = TRUE)


#
# Ptpn1 Expression
#

# Measuring Ptpn1 expression in different cell types, across treatment groups
FeaturePlot(so, features = "Ptpn1", split.by = "orig.ident", label = TRUE)


# Recreating the violin plot of PTP1B expression in AD Cells (Fig.3B)

# Creating a subset of AD cells
so_ad <- subset(so, subset = orig.ident == "AD")

# Assigning the manual cell type identities
Idents(so_ad) <- "cell_type"

# List of cell types in the same order as the paper
cell_order <- c("Choroid plexus cell", "Neuron", "Microglia", "Oligodendrocyte",
                "Astrocyte", "T cell", "OPC", "Macrophage", "Endothelial cell")
# Removing any other cell types
so_ad <- subset(so_ad, cell_type %in% cell_order)

# Sorting the cell types in order
so_ad$cell_type <- factor(so_ad$cell_type, levels = cell_order)

# Creating the violin plot
VlnPlot(so_ad, features = "Ptpn1", pt.size = 0) +
  labs(x = "Cell Types", y = "Ptpn1 Expression") +
  theme(legend.position = "none") +
  stat_summary(fun = median, geom = "point")


#
# Microglia Analysis
#

# Creating a subset containing only microglia
microglia <- subset(so, idents = "Microglia")


# Redoing the analysis on this group
microglia <- NormalizeData(microglia)

microglia <- FindVariableFeatures(microglia)
top10_microglia <- head(VariableFeatures(microglia), 10)
plot1 <- VariableFeaturePlot(microglia)
plot2 <- LabelPoints(plot = plot1, points = top10_microglia, repel = TRUE)
plot2

microglia <- ScaleData(microglia)

microglia <- RunPCA(microglia)

ElbowPlot(microglia)

microglia <- FindNeighbors(microglia, dims = 1:10)

# Lower resolution needed this time
microglia <- FindClusters(microglia, resolution = 0.2)

microglia <- RunUMAP(microglia, dims = 1:10)

DimPlot(microglia, reduction = "umap", label = TRUE)
DimPlot(microglia, reduction = "umap", label = TRUE, split.by = "orig.ident")

# Removing the incorrect clusters
microglia <- subset(microglia, idents = c("0", "1", "2"))
# Removing massive outlier cells
microglia <- subset(microglia, subset = umap_1 > -10)

# Re-plotting
DimPlot(microglia, reduction = "umap", label = TRUE)


# Identifying the microglia subclusters with marker genes

# Homeostatic microglia: Fcrls, Tmem119, P2ry12
FeaturePlot(microglia, features = "Fcrls", label = TRUE)
FeaturePlot(microglia, features = "Tmem119", label = TRUE)
FeaturePlot(microglia, features = "P2ry12", label = TRUE)
# Cluster 1

# Disease-associated microglia: Axl, Ctsl, Trem2
FeaturePlot(microglia, features = "Axl", label = TRUE)
FeaturePlot(microglia, features = "Ctsl", label = TRUE)
FeaturePlot(microglia, features = "Trem2", label = TRUE)
# Cluster 2

# Interferon-responsive microglia: Oasl2, Stat1, and Irf7
FeaturePlot(microglia, features = "Oasl2", label = TRUE)
FeaturePlot(microglia, features = "Stat1", label = TRUE)
FeaturePlot(microglia, features = "Irf7", label = TRUE)
# Cluster 0

# Cluster annotations based on these results
microglia_annotations <- c(
  "0" = "IFN",
  "1" = "Homeostatic",
  "2" = "DAM"
)
microglia@meta.data$cell_type <- microglia_annotations[as.character(microglia$seurat_clusters)]
Idents(microglia) <- "cell_type"
DimPlot(microglia, reduction = "umap", label = TRUE)
