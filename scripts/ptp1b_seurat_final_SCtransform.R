#
# Ptp1b seurat with SCTransform (and after CellBender)
#

#
# Analyzing data with Seurat
# All four samples
#

#
# Annotating the cell clusters
#

# Setting the working directory to CellBender output folder
setwd("~/ptp1b_ad/results/cellbender_output")

# Loading packages
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(harmony)

# Read in the raw matrices
KO1 <- Read10X_h5("KO_SRR36691286_cellbender_filtered_seurat.h5")
KO2 <- Read10X_h5("KO_SRR36691287_cellbender_filtered_seurat.h5")
AD1 <- Read10X_h5("AD_SRR36691288_cellbender_filtered_seurat.h5")
AD2 <- Read10X_h5("AD_SRR36691289_cellbender_filtered_seurat.h5")

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

# Filtering the low quality cells (filters from paper)
so <- subset(so, subset = nFeature_RNA > 500 & 
               nFeature_RNA < 6000 & 
               nCount_RNA > 1000 &
               nCount_RNA < 30000 & 
               percent.mt < 10)


# SCTransform, PCA, Harmony
so <- SCTransform(so, vars.to.regress = "percent.mt")
# so <- SCTransform(so, vars.to.regress = "percent.mt", variable.features.n = 2000, conserve.memory = TRUE)
so <- RunPCA(so)
so <- RunHarmony(so, group.by.vars = "orig.ident", assay.use = "SCT")


# Visualizing the PCs with an elbow plot
ElbowPlot(so)

# Clustering the cells
so <- FindNeighbors(so, dims = 1:30, reduction = "harmony")
so <- FindClusters(so, resolution = 0.15)

# Condensing all the PCs down to 2 dimensions
so <- RunUMAP(so, dims = 1:30, reduction = "harmony")

# Visualizing the clusters with a dimensional reduction plot
DimPlot(so, reduction = "umap", label = TRUE)

# Visualizing expression of one gene with a feature plot
FeaturePlot(so, features = "Ptpn1")


# SingleR #

# Identifying the marker genes (must join layers first)
so <- PrepSCTFindMarkers(so)
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
# Clusters: 3

# Choroid plexus cell: Ttr, Ecrg4, Enpp2
FeaturePlot(so, features = "Ttr", label = TRUE)
FeaturePlot(so, features = "Ecrg4", label = TRUE)
FeaturePlot(so, features = "Enpp2", label = TRUE)
# Clusters: 2, 4, 9, 12, 13, 15

# Endothelial cell: Flt1, Slco1a4, Mecom
FeaturePlot(so, features = "Flt1", label = TRUE)
FeaturePlot(so, features = "Slco1a4", label = TRUE)
FeaturePlot(so, features = "Mecom", label = TRUE)
# Clusters: 14

# Macrophage: F13a1, H2-Aa, Lyve1
FeaturePlot(so, features = "F13a1", label = TRUE)
FeaturePlot(so, features = "H2-Aa", label = TRUE)
FeaturePlot(so, features = "Lyve1", label = TRUE)
# Clusters: 11

# Microglia: P2ry12, Hexb, Trem2
FeaturePlot(so, features = "P2ry12", label = TRUE)
FeaturePlot(so, features = "Hexb", label = TRUE)
FeaturePlot(so, features = "Trem2", label = TRUE)
# Clusters: 0, 7

# Neuron: Syt1, Snap25, Nrn1
FeaturePlot(so, features = "Syt1", label = TRUE)
FeaturePlot(so, features = "Snap25", label = TRUE)
FeaturePlot(so, features = "Nrn1", label = TRUE)
# Clusters: 10

# Oligodendrocyte: Ermn, Mog, Aspa
FeaturePlot(so, features = "Ermn", label = TRUE)
FeaturePlot(so, features = "Mog", label = TRUE)
FeaturePlot(so, features = "Aspa", label = TRUE)
# Clusters: 1, 5, 16

# OPC: Pdgfra, Cacng4, Vcan
FeaturePlot(so, features = "Pdgfra", label = TRUE)
FeaturePlot(so, features = "Cacng4", label = TRUE)
FeaturePlot(so, features = "Vcan", label = TRUE)
# Clusters: 8

# T cell: Cd3d, Cd3e, Cd3g
FeaturePlot(so, features = "Cd3d", label = TRUE)
FeaturePlot(so, features = "Cd3e", label = TRUE)
FeaturePlot(so, features = "Cd3g", label = TRUE)
# Clusters: 6


# Cluster annotations based on these results
cluster_annotations <- c(
  "0" = "Microglia",
  "1" = "Oligodendrocyte",
  "2" = "Choroid plexus cell",
  "3" = "Astrocyte",
  "4" = "Choroid plexus cell",
  "5" = "Oligodendrocyte",
  "6" = "T cell",
  "7" = "Microglia",
  "8" = "OPC",
  "9" = "Choroid plexus cell",
  "10" = "Neuron",
  "11" = "Macrophage",
  "12" = "Choroid plexus cell",
  "13" = "Choroid plexus cell",
  "14" = "Endothelial cell",
  "15" = "Choroid plexus cell",
  "16" = "Oligodendrocyte",
  "17" = "Unknown",
  "18" = "Unknown"
)
so@meta.data$cell_type <- cluster_annotations[as.character(so$seurat_clusters)]
Idents(so) <- "cell_type"
so <- subset(so, cell_type != "Unknown") # Removing the unknown/difficult clusters
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
microglia <- SCTransform(microglia, vars.to.regress = "percent.mt")

microglia <- RunPCA(microglia, assay = "SCT")

microglia <- RunHarmony(microglia, group.by.vars = "orig.ident", assay.use = "SCT")

ElbowPlot(microglia)

microglia <- FindNeighbors(microglia, dims = 1:20, reduction = "harmony")

# Lower resolution needed this time
microglia <- FindClusters(microglia, resolution = 0.1)

microglia <- RunUMAP(microglia, dims = 1:20, reduction = "harmony")

DimPlot(microglia, reduction = "umap", label = TRUE)
DimPlot(microglia, reduction = "umap", label = TRUE, split.by = "orig.ident")


# Investigating cluster 2
microglia <- PrepSCTFindMarkers(microglia)
cluster2.markers <- FindMarkers(microglia, ident.1 = "2")
head(cluster2.markers, 10)




# Removing cluster 2 (choroid plexus contamination?)
microglia <- subset(microglia, idents = c("0", "1", "3"))
# Removing massive outlier cells?
#microglia <- subset(microglia, subset = umap_1 > -10)

# Rerunning after subsetting
microglia <- SCTransform(microglia, vars.to.regress = "percent.mt")
microglia <- RunPCA(microglia, assay = "SCT")
microglia <- RunHarmony(microglia, group.by.vars = "orig.ident", assay.use = "SCT")
microglia <- FindNeighbors(microglia, dims = 1:20, reduction = "harmony")
microglia <- FindClusters(microglia, resolution = 0.1)
microglia <- RunUMAP(microglia, dims = 1:20, reduction = "harmony")

# Re-plotting
DimPlot(microglia, reduction = "umap", label = TRUE)




# Identifying the microglia subclusters with marker genes

# Homeostatic microglia: Fcrls, Tmem119, P2ry12
FeaturePlot(microglia, features = "Fcrls", label = TRUE)
FeaturePlot(microglia, features = "Tmem119", label = TRUE)
FeaturePlot(microglia, features = "P2ry12", label = TRUE)
# Cluster 0

# Disease-associated microglia: Axl, Ctsl, Trem2
FeaturePlot(microglia, features = "Axl", label = TRUE)
FeaturePlot(microglia, features = "Ctsl", label = TRUE)
FeaturePlot(microglia, features = "Trem2", label = TRUE)
# Cluster 1

# Interferon-responsive microglia: Oasl2, Stat1, and Irf7
FeaturePlot(microglia, features = "Oasl2", label = TRUE)
FeaturePlot(microglia, features = "Stat1", label = TRUE)
FeaturePlot(microglia, features = "Irf7", label = TRUE)
# Cluster 2

# Cluster annotations based on these results
microglia_annotations <- c(
  "0" = "Homeostatic",
  "1" = "DAM",
  "2" = "IFN"
)
microglia@meta.data$cell_type <- microglia_annotations[as.character(microglia$seurat_clusters)]
Idents(microglia) <- "cell_type"
DimPlot(microglia, reduction = "umap", label = TRUE)
