#
# Analyzing data with Seurat
# All four samples
#

# Setting the working directory to CellRanger output folder
setwd("~/ptp1b_ad/results/cellranger_output")

# Loading packages
library(Seurat)
library(Matrix)
library(dplyr)

# Read in the raw matrices
KO1 <- Read10X("KO_SRR36691286/filtered_feature_bc_matrix/")
KO2 <- Read10X("KO_SRR36691287/filtered_feature_bc_matrix/")
AD1 <- Read10X("AD_SRR36691288/filtered_feature_bc_matrix/")
AD2 <- Read10X("AD_SRR36691289/filtered_feature_bc_matrix/")

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
# Fig.S4B

# Astrocyte: Aqp4, Acsbg1, Slc1a2
FeaturePlot(so, features = "Aqp4") # Cluster 5, 21, 22
FeaturePlot(so, features = "Acsbg1") # Cluster 5, 21, 22
FeaturePlot(so, features = "Slc1a2") # Cluster 5, 21, 22

# Choroid plexus cell: Ttr, Ecrg4, Enpp2
FeaturePlot(so, features = "Ttr", label = TRUE) # Cluster 0, 18, 14+15?
FeaturePlot(so, features = "Ecrg4", label = TRUE) # Cluster 0, 18, 14+15?
FeaturePlot(so, features = "Enpp2", label = TRUE) # Cluster 0, 18, 14+15?

# Endothelial cell: Flt1, Slco1a4, Mecom
FeaturePlot(so, features = "Flt1", label = TRUE) # Cluster 13
FeaturePlot(so, features = "Slco1a4", label = TRUE) # Cluster 13
FeaturePlot(so, features = "Mecom", label = TRUE) # Cluster 13

# Macrophage: F13a1, H2-Aa, Lyve1
FeaturePlot(so, features = "F13a1", label = TRUE) # Cluster 11, 20?
FeaturePlot(so, features = "H2-Aa", label = TRUE) # Cluster 11
FeaturePlot(so, features = "Lyve1", label = TRUE) # Cluster 11 (only some)

# Microglia: P2ry12, Hexb, Trem2
FeaturePlot(so, features = "P2ry12", label = TRUE) # Cluster 1, 4, 8, 14, 12? 22?
FeaturePlot(so, features = "Hexb", label = TRUE) # Cluster 1, 4, 8, 14, 22? 12? 
FeaturePlot(so, features = "Trem2", label = TRUE) # Cluster 1, 4, 8, 14, 12? 22? 11?

# Neuron: Syt1, Snap25, Nrn1
FeaturePlot(so, features = "Syt1", label = TRUE) # Cluster 10
FeaturePlot(so, features = "Snap25", label = TRUE) # Cluster 10, 9?
FeaturePlot(so, features = "Nrn1", label = TRUE) # Cluster 10 (only some)

# Oligodendrocyte: Ermn, Mog, Aspa
FeaturePlot(so, features = "Ermn", label = TRUE) # Cluster 2, 3, 6, 12? 15? 21?
FeaturePlot(so, features = "Mog", label = TRUE) # Cluster 2, 3, 6, 12? 15? 21?
FeaturePlot(so, features = "Aspa", label = TRUE) # Cluster 2, 3, 6, 12? 15? 21?

# OPC: Pdgfra, Cacng4, Vcan
FeaturePlot(so, features = "Pdgfra", label = TRUE) # Cluster 9
FeaturePlot(so, features = "Cacng4", label = TRUE) # Cluster 9
FeaturePlot(so, features = "Vcan", label = TRUE) # Cluster 9

# T cell: Cd3d, Cd3e, Cd3g
FeaturePlot(so, features = "Cd3d", label = TRUE) # Cluster 7, 23? 24? 25?
FeaturePlot(so, features = "Cd3e", label = TRUE) # Cluster 7, 23? 24? 25?
FeaturePlot(so, features = "Cd3g", label = TRUE) # Cluster 7, 23? 24? 25?

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
  "12" = "Unknown",
  "13" = "Endothelial cell",
  "14" = "Microglia",
  "15" = "Unknown",
  "16" = "Unknown",
  "17" = "Unknown",
  "18" = "Choroid plexus cell",
  "19" = "Unknown",
  "20" = "Unknown",
  "21" = "Unknown",
  "22" = "Unknown",
  "23" = "T cell",
  "24" = "T cell",
  "25" = "T cell",
  "26" = "Unknown"
)
so@meta.data$cell_type <- cluster_annotations[as.character(so$seurat_clusters)]
Idents(so) <- "cell_type"
DimPlot(so, reduction = "umap", label = TRUE)

# Filter out the Unknown cells
so <- subset(so, cell_type != "Unknown")
# Then re-plot
DimPlot(so, reduction = "umap", label = TRUE)







# Trying to recreate violin plot of PTP1B expression in AD Cells
# (Fig. 3B)
so_ad <- subset(so, subset = orig.ident == "AD")

Idents(so_ad) <- "SingleR.labels"

library(ggplot2)
VlnPlot(so_ad, features = "Ptpn1", pt.size = 0) +
  labs(x = "Cell Types", y = "Ptpn1 Expression") +
  theme(legend.position = "none")




