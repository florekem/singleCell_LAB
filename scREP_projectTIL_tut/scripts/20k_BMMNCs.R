# https://www.10xgenomics.com/datasets/20-k-bone-marrow-mononuclear-cells-bmmn-cs-5-ht-v-2-0-2-high-6-1-0 #nolint

library(Seurat)
library(ggplot2)
library(patchwork)
library(ProjecTILs)
library(scRepertoire)
library(tidyverse)

getwd()

# --- seurat sparse matrix ---------------------------------------------
sparse_matrix <- Read10X(
  "singleCell_LAB/scREP_projectTIL_tut/data/20k_BMMNCs/gex/"
)
names(sparse_matrix)
# [1] "Gene Expression"  "Antibody Capture"
gex <- sparse_matrix$`Gene Expression`
dim(gex)
# [1] 36601 17645
adt <- sparse_matrix$`Antibody Capture`
dim(adt)
# [1]   137 17645 [TotalSeq C Human Universal antibody cocktail]

# --- create seuratObject -----------------------------------------------------
seu <- CreateSeuratObject(
  gex,
  project = "20k_BMMNCs"
)
seu
# An object of class Seurat
# 36601 features across 17645 samples within 1 assay
# Active assay: RNA (36601 features, 0 variable features)
#  1 layer present: counts

# --- filter seuratObject ----------------------------------------------------

# --- scRep prepare T data -----------------------------------------------------

# --- scRep prepare B data -----------------------------------------------------
