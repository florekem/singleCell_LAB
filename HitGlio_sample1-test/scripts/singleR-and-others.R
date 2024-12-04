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
library(Azimuth)
library(SeuratDisk)
library(SeuratData)

?RunAzimuth

seu_singlet_azi <- RunAzimuth(seu_singlet, reference = "pbmcref")
seu_singlet_azi_cortex <- RunAzimuth(seu_singlet, reference = "humancortexref")
colnames(seu_singlet_azi@meta.data)
colnames(seu_singlet_azi_cortex@meta.data)

seu_singlet_azi <- NormalizeData(seu_singlet_azi)
Idents(seu_singlet_azi) <- "predicted.celltype.l2"
seu_singlet_azi_cortex <- NormalizeData(seu_singlet_azi_cortex)
Idents(seu_singlet_azi_cortex) <- "predicted.celltype.l2"

seu_DimPlot(
  seu_singlet_azi,
  reduction = "umap.sct",
  group.by = c("predicted.celltype.l2"),
  show = FALSE,
  save_path = "plots/scv-500-doublets/PBMCref",
  ggheight = NA,
  ggwidth = NA
)

# p1 <- DimPlot(
#   seu_singlet_azi,
#   # reduction = "umap.totalvi.sct",
#   reduction = "umap.sct",
#   pt.size = 0.5,
#   group.by = c("predicted.celltype.l2"),
#   # group.by = c("MULTI_ID", "predicted.celltype.l2"),
#   # group.by = c("seurat_clusters", "predicted.subclass"),
#   label = TRUE,
#   repel = TRUE,
#   # cols = c("CD4 TCM" = "red", "CD4 TEM" = "blue"),
#   # cols = c("CD4 TCM" = "red", "CD4 TEM" = "blue"),
#   # cols = c("CD14 Mono" = "red", "CD16 Mono" = "blue")
# )
p1
# & NoLegend()
#  & NoAxes()
table(seu_singlet_azi@meta.data$predicted.celltype.l2)

p2 <- FeaturePlot(
  seu_singlet,
  # features = c("GZMB", "PRF1"),
  # features = c("PRF1"),
  # features = c("CD8A", "CD8B"),
  features = c("CD8A"),
  # features = c("CD4", "CD8A", "CCR7", "CD62L"),
  reduction = "umap"
) & NoLegend()
VlnPlot(seurat_obj, features = c("CD4", "CD8A", "CCR7"), group.by = "clusters")

p1 | p2
p2

?DimPlot

# SaveH5Seurat(seu_singlet, "RDS/seu_singlet.h5Seurat")

## subset cluster 5 and 19 (brain cells)
Idents(seu_singlet) <- "seurat_clusters"
seu_singlet_subs <- subset(seu_singlet, idents = c("5", "19"))
seu_singlet_subs <- RunAzimuth(seu_singlet_subs, reference = "humancortexref")
seu_singlet_subs
seu_singlet_subs <- NormalizeData(seu_singlet_subs)

seu_singlet
colnames(seu_singlet@meta.data)
p1 <- DimPlot(
  seu_singlet,
  reduction = "umap",
  pt.size = 0.5,
  # group.by = c("predicted.celltype.l2"),
  # group.by = c("seurat_clusters", "predicted.celltype.l2"),
  # group.by = c("seurat_clusters", "predicted.subclass"),
  group.by = c("predicted.id"),
  label = TRUE,
  repel = TRUE
)
# & NoLegend()
#  & NoAxes()
p1



# project TILs label transfer
library(ProjecTILs)

ref_cd8 <- load.reference.map(
  "/mnt/sda4/singleCell_LAB/scREP_projectTIL_tut/data/refs/pTILs_hsa_cd8t.rds"
)

DefaultAssay(seu_singlet) <- "SCT"
seu_singlet_cd8 <- ProjecTILs.classifier(
  query = seu_singlet,
  ref = ref_cd8
)
# seu_singlet_pTILs_cd8 <- NormalizeData(seu_singlet_pTILs)

colnames(seu_singlet@meta.data)
colnames(seu_singlet_cd8@meta.data)

seu_DimPlot(
  seu_singlet_pTILs_cd4,
  reduction = "umap.sct",
  group.by = c("functional.cluster"),
  show = FALSE,
  save_path = "plots/scv-500-doublets/TILs_cd4",
  ggheight = NA,
  ggwidth = NA
)
# cols = c("CD4 TCM" = "red", "CD4 TEM" = "blue"),
# cols = c("CD8 TCM" = "red", "CD8 TEM" = "blue")
# & NoLegend()

ref_cd4 <- load.reference.map(
  "/mnt/sda4/singleCell_LAB/scREP_projectTIL_tut/data/refs/pTILs_hsa_cd4t.rds"
)
seu_singlet_pTILs_cd4 <- ProjecTILs.classifier(
  query = seu_singlet,
  ref = ref_cd4
)
# seu_singlet_pTILs_cd4 <- NormalizeData(seu_singlet_pTILs_cd4)


## DC
ref_dc <- load.reference.map(
  "/mnt/sda4/singleCell_LAB/scREP_projectTIL_tut/data/refs/pTILS_DC_human_ref_v1.rds"
)
seu_singlet_pTILs_DC <- ProjecTILs.classifier(
  query = seu_singlet,
  ref = ref_dc
)
seu_singlet_pTILs_DC <- NormalizeData(seu_singlet_pTILs_DC)

p3 <- DimPlot(
  seu_singlet_pTILs_DC,
  reduction = "umap",
  pt.size = 0.5,
  group.by = c("functional.cluster"),
  # group.by = c("seurat_clusters", "predicted.celltype.l2"),
  # group.by = c("seurat_clusters", "predicted.subclass"),
  label = TRUE,
  repel = TRUE,
  # cols = c("CD4 TCM" = "red", "CD4 TEM" = "blue"),
  # cols = c("CD8 TCM" = "red", "CD8 TEM" = "blue")
)
p3

## subset cluster with tumor CD8 and transfer labels from TILs
DimPlot(
  seu_singlet_pTILs_DC,
  reduction = "umap",
  pt.size = 0.5,
  group.by = c("seurat_clusters"),
  # group.by = c("seurat_clusters", "predicted.celltype.l2"),
  # group.by = c("seurat_clusters", "predicted.subclass"),
  label = TRUE,
  repel = TRUE,
  # cols = c("CD4 TCM" = "red", "CD4 TEM" = "blue"),
  # cols = c("CD8 TCM" = "red", "CD8 TEM" = "blue")
)
seu_singlet_subs_tumor_cd8 <- subset(seu_singlet, idents = c("5", "19"))

plot <- FeaturePlot(
  seu_singlet_pTILs,
  # features = c("CD69", "CD103", "CD49a"),
  # features = c("PD1", "CTLA4", "TIM3", "LAG3", "TIGIT"),
  # features = c("CD28", "CD27", "IL7R"),
  # features = c("EOMES", "TOX", "BCL6", "CXCR5"),
  # features = c("Ki67", "PCNA", "MCM6"),
  # features = c("TMEM119", "P2RY12", "CSF1R", "CD11b", "CD45"),
  # features = c("CD11c", "CD80", "CD86", "MHCII", "CD40", "CLEC9A", "ITGAX"),
  # features = c("CLEC9A", "CD40", "MHCII", "ITGAX", "CD11c", "CD80", "CD83"),
  # features = c("CD14", "CD68", "CD11b", "CSF1R", "CD163", "MRC1", "MARCO"),
  # features = c("CX3CR1", "TREM1", "TREM2", "P2RY12", "CSF1R", "TMEM119", "ITGAX", "CD83", "CD14"),
  # features = c("CD8A", "CD4", "TCRB", "CD3D", "CD3E", "CD86", "ITGAX", "CD68", "CD163", "CD14", "CD11b", "ITGAE"),
  features = "IL7R",
  reduction = "umap",
  label = TRUE
)
ggsave("plots/IL7R.png", plot)
plot

macrophage_markers <- c(
  "HLA-DRA", "HLA-DRB1", "CD80", "CD86", "CD163", "CD14",
  "IRF5", "NOS2", "CXCL10", "CCL2", "TGFB1", "IL1B", "IL6"
)
cd8_t_cell_markers <- c(
  "CD8B", "GZMK", "GZMA", "PRF1", "IFNG", "TNF", "CD28",
  "CD27", "IL7R", "TIGIT", "LAG3"
)
