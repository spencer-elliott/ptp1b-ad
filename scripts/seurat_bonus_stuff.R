
# Endothelial Cell Analysis ----

# Loading the Seurat object
so <- readRDS("~/ptp1b_ad/results/seurat_objects/ptp1b_so_final.rds")

# Creating a subset containing only endothelial cells
endothelial <- subset(so, idents = "Endothelial cell")

# Fcgrt expression in endothelial cells
FeaturePlot(endothelial, features = "Fcgrt", split.by = "orig.ident", label = TRUE)

# Investigating differentially expressed genes
endothelial <- PrepSCTFindMarkers(endothelial)
DefaultAssay(endothelial) <- "RNA"
endothelial <- NormalizeData(endothelial)
endothelial.markers <- FindMarkers(endothelial, ident.1 = "AD", ident.2 = "KO", group.by = "orig.ident")
# Top 20 differentially expressed genes
head(endothelial.markers, 20)
# Fcgrt expression
endothelial.markers["Fcgrt", ]

saveRDS(endothelial, "~/ptp1b_ad/results/seurat_objects/ptp1b_endothelial_final.rds")
# To load: endothelial <- readRDS("~/ptp1b_ad/results/seurat_objects/ptp1b_endothelial_final.rds")


# B Cell Analysis ----

# Subset of cells labeled as B cells by SingleR
b_cell <- subset(so, subset = SingleR.labels == "B cells")
# And macrophage by manual annotations
b_cell <- subset(b_cell, idents = "Macrophage")

# Cd19 expression in B cells
FeaturePlot(b_cell, features = "Cd19", split.by = "orig.ident", label = TRUE)
# Cd79a expression in B cells
FeaturePlot(b_cell, features = "Cd79a", split.by = "orig.ident", label = TRUE)
# Prdm1 expression in B cells
FeaturePlot(b_cell, features = "Prdm1", split.by = "orig.ident", label = TRUE)

