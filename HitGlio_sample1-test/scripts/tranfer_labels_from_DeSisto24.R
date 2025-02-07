# this data comes from


# seu_singlet from 01_demux_filter.R

# 1. Paths ------------------------------------------------------
sample_path <- file.path("/mnt/sda4/PediatricNeuroOncologyCellAtlas/DeSisto24")
print(sample_path)

seu_singlet <- readRDS("RDS/seu_singlet_clustered_totalvi-28nov2024.rds")

# 2. Create seurat object ---------------------------------------
# test <- read_tsv(file.path(sample_path, "tumor", "exprMatrix.tsv"))
test <- read_tsv(file.path(sample_path, "immune", "exprMatrix.tsv"))
head(test, n = 1)
gene_names <- test$gene
test <- test[, -1]
test_mtx <- as.matrix(test)
rownames(test)
rownames(test_mtx) <- gene_names
rownames(test_mtx)
colnames(test_mtx)

seu_test <- CreateSeuratObject(
  counts = test_mtx,
  project = "thal_single_cell",
  min.cells = 3,
  min.features = 200
)
seu_test
head(seu_test@meta.data)
unique(seu_test@meta.data$orig.ident)
table(seu_test@meta.data$orig.ident)

seu <- seu_test

# annot <- read.csv(file.path(sample_path, "tumor", "meta.tsv"))
annot <- read.csv(file.path(sample_path, "immune", "meta.tsv"), sep = "\t")

head(annot)
rownames(annot) <- annot$Cell

## 2.1. filtering by annotation file --------
seu <- subset(
  seu,
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
Idents(seu) <- "cell.types"
unique(Idents(seu))
table(Idents(seu))

head(seu@meta.data)
Features(seu)
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
  reduction = "pca", dims = 1:10
)
seu <- FindClusters(
  seu,
  resolution = 0.5, verbose = FALSE,
  algorithm = 4
)
seu <- RunUMAP(
  seu,
  reduction = "pca", dims = 1:10
)
DimPlot(
  seu,
  group.by = "sample", reduction = "umap",
  label = TRUE
)



colnames(seu@meta.data)
## 3.2 add sample name metadata (for scvi integration) ----------------
seu$data_name <- "DeSisto24"
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
  dims = 1:10,
  reference.reduction = "pca"
)
colnames(integrated@meta.data)
predictions <- TransferData(
  anchorset = transf.anchors,
  refdata = integrated$cell.types,
  dims = 1:10
)

seu_singlet_myelo <- AddMetaData(
  seu_singlet_myelo,
  metadata = predictions
)

DefaultAssay(seu_singlet) <- "RNA"
Idents(seu_singlet_myelo) <- "predicted.id"
seu_singlet <- NormalizeData(seu_singlet)
seu_singlet <- FindVariableFeatures(seu_singlet)
seu_singlet <- ScaleData(seu_singlet)
seu_singlet <- RunPCA(seu_singlet)

DimPlot(
  seu_singlet_myelo,
  # group.by = c("my_id", "predicted.id"),
  group.by = c("predicted.id"),
  reduction = "wnn.umap",
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
