# this data comes from


# seu_singlet from 01_demux_filter.R

seu_singlet <- readRDS("/mnt/sda4/singleCell_LAB/HitGlio_sample1-test/RDS/06-dec-24_seu_singlet_final_doublet-500_wo-cd8-myelo_sct_totalvi.rds")

seu_singlet <- seu_normalize_var_scale_pca(
  seu_singlet,
  normalize = "log",
  cluster.name = "cluster.pca",
  run_pca = TRUE,
  pca.reduction.name = "pca",
  umap.reduction.name = "umap",
  dims = 1:20,
  k.param = 20,
  algorithm = 4,
  resolution = 0.5,
  group.singletons = FALSE
)


# 1. Paths ------------------------------------------------------
sample_path <- file.path("/mnt/sda4/brain_immune_atlas/humanRecurrentGBMfullAggregate/filtered_matrix/")
print(sample_path)

# 2. Create seurat object ---------------------------------------
# test <- read_tsv(file.path(sample_path, "tumor", "exprMatrix.tsv"))

sparse_matrix <- Seurat::Read10X(
  data.dir = sample_path
)

seu_test <- CreateSeuratObject(
  counts = sparse_matrix,
  project = "recurrentGBM",
  min.cells = 3,
  min.features = 200
)

seu_test
head(seu_test@meta.data)
unique(seu_test@meta.data$orig.ident)
table(seu_test@meta.data$orig.ident)

seu <- seu_test

# annot <- read.csv(file.path(sample_path, "tumor", "meta.tsv"))
annot <- read.csv(file.path(sample_path, "/../", "annot_Human_R_GBM_Full.csv"))

head(annot)
rownames(annot) <- annot$cell

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
Idents(seu) <- "cluster"
unique(Idents(seu))
table(Idents(seu))

head(seu@meta.data)

# 3.1 process new dataset (prepare to integrate) ------------------
# Integration with log-normalized data
# it does not work with SCT data?
# seu <- NormalizeData(seu)
# seu <- FindVariableFeatures(seu)
# seu <- ScaleData(seu)
# seu <- RunPCA(seu)
# ElbowPlot(seu, ndims = 50)
# #
# seu <- FindNeighbors(
#   seu,
#   reduction = "pca", dims = 1:10
# )
# seu <- FindClusters(
#   seu,
#   resolution = 0.5, verbose = FALSE,
#   algorithm = 4
# )
# seu <- RunUMAP(
#   seu,
#   reduction = "pca", dims = 1:10
# )
seu <- seu_normalize_var_scale_pca(
  seu,
  normalize = "log",
  cluster.name = "cluster.pca",
  run_pca = TRUE,
  pca.reduction.name = "pca",
  umap.reduction.name = "umap",
  dims = 1:20,
  k.param = 20,
  algorithm = 4,
  resolution = 0.5,
  group.singletons = FALSE
)

DimPlot(
  seu,
  group.by = "cluster", reduction = "umap.pca.sct",
  label = TRUE
)

## 3.2 add sample name metadata (for scvi integration) ----------------
seu$data_name <- "DeSisto24"
seu_singlet$data_name <- "sample1"
head(seu@meta.data)

# 4. Merge both datasets --------------------------
DefaultAssay(seu) <- "RNA"
DefaultAssay(seu_singlet) <- "RNA"

merged <- merge(seu_singlet, seu)
merged
merged <- seu_normalize_var_scale_pca(
  merged,
  normalize = "log",
  cluster.name = "cluster.pca",
  run_pca = TRUE,
  pca.reduction.name = "pca",
  umap.reduction.name = "umap",
  dims = 1:20,
  k.param = 20,
  algorithm = 1,
  resolution = 0.5,
  group.singletons = FALSE
)
# merged <- NormalizeData(merged)
# merged <- FindVariableFeatures(merged)
# merged <- ScaleData(merged)
# merged <- RunPCA(merged)

# 5. Integrate ----------------------
integrated <- IntegrateLayers(
  merged,
  # method = CCAIntegration,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  assay = "RNA",
  new.reduction = "integrated.harmony"
)

## 5.1 Inspect intergration ------------------
integrated <- FindNeighbors(
  integrated,
  reduction = "integrated.harmony",
  dims = 1:10
)
integrated <- FindClusters(
  integrated,
  resolution = 0.3,
  algorithm = 1
)
integrated <- RunUMAP(
  integrated,
  reduction = "integrated.harmony", dims = 1:10,
  reduction.name = "umap.integrated"
)
colnames(integrated@meta.data)
DimPlot(
  integrated,
  group.by = c("MULTI_ID"),
  reduction = "umap.integrated",
  label = TRUE
)
# 6. Transfer labels -------------------
transf.anchors <- FindTransferAnchors(
  reference = integrated,
  query = seu_singlet,
  dims = 1:10,
  reference.reduction = "pca"
)
colnames(integrated@meta.data)

predictions <- TransferData(
  anchorset = transf.anchors,
  refdata = integrated$cluster,
  dims = 1:10
)

seu_singlet <- AddMetaData(
  seu_singlet,
  metadata = predictions
)

colnames(seu_singlet@meta.data)
seu_singlet
seu_DimPlot(
  seu_singlet,
  reduction = "umap.totalvi.sct",
  # group.by = c("cluster.totalvi", "predicted.id"),
  group.by = c("MULTI_ID", "cluster.totalvi"),
  show = TRUE,
  # save_path = "plots/clusters",
  ggheight = NA,
  ggwidth = NA,
  label.size = 4
)
















DefaultAssay(seu_singlet) <- "RNA"
Idents(seu_singlet) <- "predicted.id"
seu_singlet <- NormalizeData(seu_singlet)
seu_singlet <- FindVariableFeatures(seu_singlet)
seu_singlet <- ScaleData(seu_singlet)
seu_singlet <- RunPCA(seu_singlet)

DimPlot(
  seu_singlet,
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
