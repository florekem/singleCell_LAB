# this data comes from


# seu_singlet from 01_demux_filter.R

# 1. Paths ------------------------------------------------------
sample_path <- file.path("/mnt/sda4/brain_immune_atlas/human-gbm-citeseq-TAM/")
print(sample_path)

seu_singlet <- readRDS("RDS/seu_singlet_clustered_totalvi-28nov2024.rds")

# 2. Create seurat object ---------------------------------------
seu_matrix <- Read10X(file.path(sample_path, filtered_feature_bc_matrix))

seu_test <- CreateSeuratObject(
  counts = seu_matrix$`Gene Expression`,
  project = "thal_single_cell",
)
seu_test
head(seu_test@meta.data)
unique(seu_test@meta.data$orig.ident)
table(seu_test@meta.data$orig.ident)

seu <- seu_test

# annot <- read.csv(file.path(sample_path, "tumor", "meta.tsv"))
annot <- read.csv(file.path(sample_path, "annot_Human_TAM_DC_Mono_citeSeq.csv"))
head(annot)
rownames(annot) <- annot$cell

## 2.1. filtering by annotation file --------
seu <- subset(
  seu_test,
  cells = rownames(annot)
)
seu

## 2.2. add metadata from anot file ---------------------
seu <- AddMetaData(
  seu,
  metadata = annot
)
seu
colnames(seu@meta.data)
Idents(seu) <- "cluster"
unique(Idents(seu))
table(Idents(seu))

head(seu@meta.data)

# 3.1 process new dataset (prepare to integrate) ------------------
# Integration with log-normalized data
# it does not work with SCT data?
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
ElbowPlot(seu, ndims = 50)

seu <- FindNeighbors(
  seu,
  reduction = "pca", dims = 1:30
)
seu <- FindClusters(
  seu,
  resolution = 0.5, verbose = FALSE,
  algorithm = 4
)
seu <- RunUMAP(
  seu,
  reduction = "pca", dims = 1:30
)
DimPlot(
  seu,
  group.by = "cluster", reduction = "umap",
  label = TRUE
)

# subset only myelo clusters from seu_singlet
seu_singlet_myelo <- subset(seu_singlet, idents = c(1, 11, 10, 18, 15, 20, 24, 29, 28, 9, 16, 27, 12))

### lets try without integration
seu_singlet_myelo <- NormalizeData(seu_singlet_myelo)
seu_singlet_myelo <- FindVariableFeatures(seu_singlet_myelo)
seu_singlet_myelo <- ScaleData(seu_singlet_myelo)
seu_singlet_myelo <- RunPCA(seu_singlet_myelo)
ElbowPlot(seu_singlet_myelo, ndims = 50)

seu_singlet_myelo <- FindNeighbors(
  seu_singlet_myelo,
  reduction = "pca", dims = 1:30
)
seu_singlet_myelo <- FindClusters(
  seu_singlet_myelo,
  resolution = 0.5, verbose = FALSE,
  algorithm = 4
)
seu_singlet_myelo <- RunUMAP(
  seu_singlet_myelo,
  reduction = "pca", dims = 1:30
)
DimPlot(
  seu_singlet_myelo,
  group.by = "seurat_clusters", reduction = "umap",
  label = TRUE
)
colnames(seu_singlet_myelo@meta.data)

anchors <- FindTransferAnchors(
  reference = seu,
  query = seu_singlet_myelo, dims = 1:30
)

query <- TransferData(
  anchorset = anchors,
  refdata = seu$cluster, dims = 1:30
)

seu_singlet_myelo$predicted.id <- query$predicted.id
seu_singlet_myelo$prediction.score <- query$prediction.score.max
seu_singlet_myelo$filtered_prediction <- ifelse(
  seu_singlet_myelo$prediction.score > 0.8,
  seu_singlet_myelo$predicted.id, "Uncertain"
)

DimPlot(
  seu_singlet_myelo,
  group.by = c("seurat_clusters", "filtered_prediction"),
  label = TRUE, repel = TRUE,
  # cols = c("Hypoxic Mo-TAM" = "red")
  cols = c("Classical monocytes" = "red")
)

### tu chciałbym spórbować przenieść powyższe (na samych mylo)
### od razu na seu_singlet

all(colnames(seu_singlet_myelo) %in% colnames(seu_singlet))

myeloid_predictions <- data.frame(
  cell = rownames(seu_singlet_myelo@meta.data),
  predicted.id = seu_singlet_myelo$filtered_prediction,
  prediction.score = seu_singlet_myelo$prediction.score
)

# Create a data frame for all cells with default NA values
all_predictions <- data.frame(
  cell = colnames(seu_singlet),
  predicted.id = NA,
  prediction.score = NA
)

# Update the predictions for the myeloid subset
all_predictions <- merge(
  all_predictions, myeloid_predictions,
  by = "cell", all.x = TRUE
)

# Add the predictions back to the full Seurat object
seu_singlet$predicted.id <- all_predictions$predicted.id
seu_singlet$prediction.score <- all_predictions$prediction.score


DimPlot(
  seu_singlet,
  reduction = "umap.totalvi.sct",
  group.by = c("cluster.totalvi", "predicted.id"),
  label = TRUE, repel = TRUE,
  cols = c("Hypoxic Mo-TAM" = "red")
  # cols = c("Classical monocytes" = "red")
  # cols = c("Transitory Mo-TAM" = "red")
  # cols = c("mg-TAM" = "red")
  # cols = c("IFN Mo-TAM" = "red")
)

FeaturePlot(
  seu_singlet,
  reduction = "umap.totalvi.sct",
  # features = c("FCN1", "VCAN", "CD14", "FCGR3A")
  # features = c("TGFBI", "CLEC12A", "FXYD5")
  # features = c("SALL1", "TMEM119", "P2RY12", "TMIGD3", "APOC2", "SCIN")
  # features = c("BNIP3", "ADAM8", "MIF", "SLC2A1", "ADM", "CSTB")
  # features = c("CX3CR1", "BIN1", "SCIN")
  # features = c("OLFML3", "MSR1", "CD163", "OLR1", "C3", "LPCAT2", "CXCL16")
  # features = c("GFAP")
)


pca_loadings <- Loadings(seu[["pca"]])
top_genes <- head(rownames(pca_loadings[order(abs(pca_loadings[, 1]), decreasing = TRUE), ]), 10)
top_genes

Idents(seu) <- "cluster"
markers <- FindMarkers(
  seu,
  ident.1 = "Hypoxic Mo-TAM",
  test_use = "MAST",
  min.pct = 0.25
)
markers[1:50, ]
check %in% rownames(markers)

common <- cluster_markers[rownames(cluster_markers) %in% rownames(markers), ]
dim(common)
common[1:50, ]

separate <- cluster_markers[!rownames(cluster_markers) %in% rownames(markers), ]
dim(separate)

cluster_markers <- FindMarkers(
  seu_singlet,
  test_use = "MAST",
  ident.1 = 1,
  min.pct = 0.25
)
rownames(cluster_markers)[1:50]

cluster_markers[1:50, ]

check <- c("BNIP3", "ADAM8", "MIF", "SLC2A1")
check %in% rownames(cluster_markers)



### check for ModuleScore
seu_module_score <- Seurat::AddModuleScore(
  seu_singlet,
  features = list(hypoxia_tam),
  name = "hypoxia_tam"
)
hypoxia_tam <- c(
  "ENO1", "SCD", "PLP2", "RAB42", "S100A10", "S100A6", "HK2", "SLC2A1",
  "CXCL8", "LGALS1", "TIMP1", "PLIN2", "CTSL", "LDHA", "NDRG1", "HILPDA",
  "ERO1A", "NUPR1", "MT2A"
)
Seurat::VlnPlot(
  seu_module_score,
  "hypoxia_tam1"
)
FeaturePlot(
  seu_module_score,
  "hypoxia_tam1",
  reduction = "umap.totalvi.sct"
)
?FeaturePlot

head(seu_module_score@meta.data)
dittoSeq::dittoBarPlot(seu_module_score, var = "filtered_prediction", group.by = "cluster.totalvi")


## 3.2 add sample name metadata (for scvi integration) ----------------
seu$data_name <- "gbm-citeseq"
seu_singlet_myelo$data_name <- "sample1"
head(seu@meta.data)

# 4. Merge both datasets --------------------------
DefaultAssay(seu) <- "RNA"
DefaultAssay(seu_singlet_myelo) <- "RNA"

merged <- merge(seu_singlet_myelo, seu)
merged
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)

# 5. Integrate ----------------------
integrated <- IntegrateLayers(
  merged,
  # method = CCAIntegration,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  assay = "RNA",
  new.reduction = "integrated.cca"
)
?IntegrateLayers

## 5.1 Inspect intergration ------------------
integrated <- FindNeighbors(
  integrated,
  reduction = "integrated.cca",
  dims = 1:10
)
integrated <- FindClusters(
  integrated,
  resolution = 0.3,
  algorithm = 1
)
integrated <- RunUMAP(
  integrated,
  reduction = "integrated.cca", dims = 1:10,
  reduction.name = "umap.integrated"
)
colnames(integrated@meta.data)
DimPlot(
  integrated,
  group.by = c("data_name"),
  reduction = "umap.integrated",
  label = TRUE
)
# 6. Transfer labels -------------------
transf.anchors <- FindTransferAnchors(
  reference = integrated,
  query = seu_singlet_myelo,
  dims = 1:5,
  reference.reduction = "pca"
)
colnames(integrated@meta.data)

predictions <- TransferData(
  anchorset = transf.anchors,
  refdata = integrated$cluster,
  dims = 1:5
)

seu_singlet_myelo <- AddMetaData(
  seu_singlet_myelo,
  metadata = predictions
)

DefaultAssay(seu_singlet_myelo) <- "RNA"
Idents(seu_singlet_myelo) <- "predicted.id"
seu_singlet_myelo <- NormalizeData(seu_singlet_myelo)
seu_singlet_myelo <- FindVariableFeatures(seu_singlet_myelo)
seu_singlet_myelo <- ScaleData(seu_singlet_myelo)
seu_singlet_myelo <- RunPCA(seu_singlet_myelo)
seu_singlet_myelo <- RunUMAP(
  seu_singlet_myelo,
  reduction = "pca", dims = 1:5,
  reduction.name = "umap"
)

DimPlot(
  seu_singlet_myelo,
  # group.by = c("my_id", "predicted.id"),
  group.by = c("predicted.id"),
  reduction = "umap",
  label = TRUE
)
table(seu_singlet@meta.data$predicted.id)

?DimPlot
seu_singlet

integrated <- RunUMAP(
  integrated,
  dims = 1:30,
  reduction = "integrated.harmony",
  return.model = TRUE
)

mapquer <- MapQuery(
  anchorset = transf.anchors,
  reference = integrated,
  query = seu_singlet,
  refdata = list(celltype = "cluster"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

p1 <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "cluster", label = TRUE, label.size = 3,
  repel = TRUE
) + NoLegend() + ggtitle("Reference annotations")

p2 <- DimPlot(
  seu_singlet,
  reduction = "umap",
  group.by = "predicted.id",
  label = TRUE,
  label.size = 3, repel = TRUE
) + NoLegend() + ggtitle("Query transferred labels")

p1 + p2

# 7. Integrate using scVI (replaces 6.) ---------------------------------------
# /saveRDS(merged, "RDS/merged-temp.rds")
merged <- readRDS("RDS/merged-temp.rds")
library(sceasy)
library(reticulate)

merged <- JoinLayers(merged)
merged

## /issues with v5
merged[["RNA3"]] <- as(object = merged[["RNA"]], Class = "Assay")
DefaultAssay(merged) <- "RNA3"
merged[["RNA"]] <- NULL
merged <- RenameAssays(object = merged, RNA3 = "RNA")

adata <- convertFormat(
  merged,
  from = "seurat",
  to = "anndata",
  main_layer = "counts",
  drop_single_values = FALSE,
  outFile = "RDS/filename.h5ad"
)
adata

## Run this in python --------------------------
# run setup_anndata, use column stim for batch
scvi$model$SCVI$setup_anndata(adata, batch_key = "MULTI_ID")
# create the model
model <- scvi$model$SCVI(adata)
# train the model
model$train()
# to specify the number of epochs when training:
# model$train(max_epochs = as.integer(400))

reticulate::use_condaenv("/mnt/sda4/conda/envs/R")

use_condaenv("R")

## And go back to R with latent.csv ---------------

latent <- read.csv("RDS/latent.csv")
latent_mtx <- as.matrix(latent)
latent_mtx
rownames(latent_mtx) <- colnames(merged)
latent_mtx <- latent_mtx[, -1] # remove 1st col
colnames(latent_mtx) <- paste0("scvi_", 1:10)
colnames(latent_mtx)

merged[["scvi"]] <- CreateDimReducObject(
  embeddings = latent_mtx,
  key = "scvi_",
  assay = DefaultAssay(merged)
)

merged <- FindNeighbors(merged, dims = 1:10, reduction = "scvi")
merged <- FindClusters(merged, resolution = 1)
merged <- RunUMAP(merged, dims = 1:10, reduction = "scvi", n.components = 2)
DimPlot(merged, reduction = "umap", pt.size = 3)
DimPlot(
  merged,
  reduction = "umap",
  pt.size = 0.5,
  group.by = "cluster",
  label = TRUE
)



# baza
seu_singlet_myelo
DefaultAssay(seu_singlet_myelo) <- "RNA"

top_genes <- Seurat::FindAllMarkers(
  seu_singlet_myelo,
  test_use = "MAST", only.pos = TRUE
) |>
  group_by(cluster) |>
  top_n(n = 5, wt = avg_log2FC)

DoHeatmap(seu_singlet_myelo, features = unique(top_genes$gene)) +
  scale_fill_viridis_c()


Idents(seu_singlet_myelo) <- "cluster.totalvi"
markers <- FindMarkers(
  seu_singlet_myelo,
  ident.1 = "10",
  test_use = "MAST",
  min.pct = 0.25
) |>
  mutate(pct1pct2 = pct.1 - pct.2) |>
  # filter(pct1pct2 > 0.5) |>
  # filter(pct1pct2 > 0.3) |>
  arrange(desc(avg_log2FC))

dim(markers)
markers[1:50, ]

"EREG" %in% rownames(markers)
"IFITM2" %in% rownames(markers)
microglia <- c("TMEM119", "P2RY12", "CX3CR1", "P2RY13", "SELPLG", "CSF1R", "TREM2", "CCR5", "MARCKS")
test_this <- microglia
test_this %in% rownames(markers)
markers[test_this, ]

FeaturePlot(
  seu_singlet_myelo,
  reduction = "umap.totalvi.sct",
  features = "TMEM119"
)

"P2YR12" %in% markers


DimPlot(
  seu_singlet_myelo,
  group.by = "cluster.totalvi", reduction = "umap.totalvi.sct",
  label = TRUE
)


markers <- FindMarkers(
  seu_singlet_myelo,
  ident.1 = "18",
  ident.2 = "1",
  test_use = "MAST",
  min.pct = 0.25
) |>
  mutate(pct1pct2 = pct.1 - pct.2) |>
  # filter(pct1pct2 > 0.5) |>
  filter(pct1pct2 > 0.3) |>
  arrange(desc(avg_log2FC))
dim(markers)
markers[1:50, ]
rownames(markers)
