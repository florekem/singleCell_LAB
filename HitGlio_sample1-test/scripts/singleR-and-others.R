library(celldex)
library(SingleR)
library(scRNAseq)
library(scuttle)
library(sceasy)

setwd("/mnt/sda4/singleCell_LAB/HitGlio_sample1-test")
seu_singlet <- readRDS("RDS/seu_singlet_after-filtering.rds")

# 1. sce ---------------------------------------------------------
# /something not working here?
searchDatasets("brain")[, c("name", "title", "version", "genome")]

sce <- fetchDataset("darmanis-brain-2015", "2023-12-21")
sce
sce <- sce[, !is.na(sce$label)] # remove unlabeled
sce@assays
sce <- logNormCounts(sce)

pred.grun <- SingleR(test = sceG, ref = sceM, labels = sceM$label, de.method = "wilcox")

sce_singlet <- as.SingleCellExperiment(seu_singlet)
sce_singlet

seu_SingleR_sce <- SingleR::SingleR(
  test = Seurat::GetAssayData(seu_singlet),
  ref = sce,
  labels = sce$label
)

## 1.1. convert sce to seurat and try to integrate with seu_singlet ----------
# https://github.com/satijalab/seurat/issues/6852
assay(sce, "counts") <- as.matrix(assay(sce, "counts"))
assay(sce, "logcounts") <- as.matrix(assay(sce, "logcounts"))

sce_to_seurat <- as.Seurat(
  sce,
  counts = "counts"
)
sce_to_seurat
colnames(sce_to_seurat@meta.data)
unique(sce_to_seurat@meta.data$cell.type)
# [1] "oligodendrocytes"  "hybrid"            "astrocytes"
# [4] "OPC"               "microglia"         "neurons"
# [7] "endothelial"       "fetal_quiescent"   "fetal_replicating"
seu <- sce_to_seurat

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
  resolution = 1.0, verbose = FALSE,
  algorithm = 4
)
seu <- RunUMAP(
  seu,
  reduction = "pca", dims = 1:30
)
DimPlot(
  seu,
  group.by = c("seurat_clusters", "cell.type"),
  reduction = "umap",
  label = TRUE
)

## 1.1 summary --------------------
# /in summary: this dataset is not worth it, but the script is here
# /for the future use

# 2. celldex ------------------------------------------------------------------
celdex <- celldex::NovershternHematopoieticData()
class(celdex)
table(celdex$label.main)

seu_SingleR <- SingleR::SingleR(
  test = Seurat::GetAssayData(seu_singlet),
  ref = celdex,
  labels = celdex$label.main
)
seu_SingleR
head(seu_SingleR)

seu_singlet$signleR_anot <- seu_SingleR$labels

SingleR::plotScoreHeatmap(seu_SingleR)

SingleR::plotDeltaDistribution(seu_SingleR)

Seurat::DimPlot(
  seu_singlet,
  group.by = "signleR_anot",
  label = TRUE
)

# 3. Azimuth --------------------------------
# not working yet!
library(Azimuth)
RunAzimuth()


library(SeuratDisk)
library(SeuratData)

SaveH5Seurat(seu_singlet, "RDS/seu_singlet.h5Seurat")
