library(Seurat)
library(ggplot2)
library(patchwork)
library(ProjecTILs)
library(scRepertoire)
library(tidyverse)

sample <- read.csv("/mnt/sda4/singleCell_LAB/scREP_projectTIL_tut/data/2k_mouse/2k_BEAM-T_Mouse_H2Kb_OT-1_5pv2_2k_BEAM-T_Mouse_H2Kb_OT-1_5pv2_vdj_t_filtered_contig_annotations.csv") # nolimes

sample_list <- list(sample)

combined_tcr <- combineTCR(
  sample_list,
  # to make things easier, w/o 'samples' argument
  # it wont add samples names to barcodes
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE
)
combined_tcr[[1]]$barcode
# "TTTGGTTCATTGTGCA-1"

clonalQuant(
  combined_tcr,
  cloneCall = "strict",
  chain = "both",
  scale = TRUE
)

clonalLength(
  combined_tcr,
  cloneCall = "aa",
  chain = "both"
)

# --- Seurat gex -------------------------------------------------
sparse_matrix <- Read10X(
  "singleCell_LAB/scREP_projectTIL_tut/data/2k_mouse/gex"
)
sparse_matrix
gex <- sparse_matrix$`Gene Expression`
adt <- sparse_matrix$`Antigen Capture`

seu <- CreateSeuratObject(
  counts = gex,
  project = "2k_mouse",
  min.cells = 5,
  min.features = 300
)
seu
# An object of class Seurat
# 13624 features across 1719 samples within 1 assay
# Active assay: RNA (13624 features, 0 variable features)
#  1 layer present: counts
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- FindNeighbors(
  seu,
  reduction = "pca", dims = 1:20
)
seu <- FindClusters(
  seu,
  resolution = 0.5, verbose = FALSE
)
seu <- RunUMAP(
  seu,
  reduction = "pca", dims = 1:20
)

# --- combibe tcr ----------------------------------------------------
seu_tcr <- combineExpression(
  combined_tcr,
  seu,
  cloneCall = "gene",
  proportion = TRUE
)
colnames(seu_tcr@meta.data)
#  [1] "orig.ident"    "nCount_RNA"       "nFeature_RNA"     "CTgene"
#  [5] "CTnt"         "CTaa"             "CTstrict"         "clonalProportion"
#  [9] "clonalFrequency"  "cloneSize"
colorblind_vector <- hcl.colors(n = 7, palette = "inferno", fixup = TRUE)

DimPlot(
  seu_tcr,
  group.by = "cloneSize"
) + scale_color_manual(values = rev(colorblind_vector[c(1, 3, 4, 5, 7)]))
