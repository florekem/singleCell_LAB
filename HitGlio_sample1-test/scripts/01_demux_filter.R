# according to: https://satijalab.org/seurat/articles/hashing_vignette

library(Seurat)
library(ggplot2)
library(Matrix)
library(tidyverse)
library(clustree)

setwd("/mnt/sda4/singleCell_LAB/HitGlio_sample1-test")

# prepare list of Seurat objects.
# seu_list <- list()
# sample_names <- c("sample1")

# sample_path <- file.path("data", "sample1")
sample_path <- file.path("/mnt/sdb1/runs/sample1_multilane_spec_index_NNN/outs/per_sample_outs/sample1/count/sample_filtered_feature_bc_matrix")
print(sample_path)

sparse_matrix <- Seurat::Read10X(data.dir = sample_path)
gex <- sparse_matrix$`Gene Expression`
adt <- sparse_matrix$`Antibody Capture`[grepl(
  "TotalSeqC", rownames(sparse_matrix$`Antibody Capture`)
), ]
hto <- sparse_matrix$`Antibody Capture`[grepl(
  "Hashtag", rownames(sparse_matrix$`Antibody Capture`)
), ]
rownames(hto)
rownames(hto) <- c("tumor", "csf")
joint_bcs <- intersect(colnames(gex), colnames(hto))
head(joint_bcs)

gex <- gex[, joint_bcs]
adt <- adt[, joint_bcs]

hto <- as.matrix(hto[, joint_bcs])

# adt <- adt[, joint_bcs]
# ncol(adt)
# hto <- hto[, joint_bcs]
# ncol(hto)

seu <- CreateSeuratObject(
  counts = Matrix::Matrix(as.matrix(gex), sparse = TRUE)
)

# Normalize RNA data with log normalization
seu <- NormalizeData(seu)
# Find and scale variable features
seu <- FindVariableFeatures(seu, selection.method = "mean.var.plot")
seu <- ScaleData(seu, features = VariableFeatures(seu))

# Add HTO data as a new assay independent from RNA
seu[["HTO"]] <- CreateAssayObject(counts = hto)
seu[["ADT"]] <- CreateAssayObject(counts = adt)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
seu <- NormalizeData(seu, assay = "HTO", normalization.method = "CLR")

seu_demux <- MULTIseqDemux(
  seu,
  assay = "HTO",
  quantile = 0.7,
  autoThresh = TRUE,
  maxiter = 5,
  qrange = seq(from = 0.1, to = 0.9, by = 0.05),
  verbose = TRUE
)
table(seu_demux$MULTI_ID)
table(seu_demux$MULTI_classification)

Idents(seu_demux) <- "MULTI_ID"
VlnPlot(seu_demux, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

RidgePlot(
  seu_demux,
  assay = "HTO",
  features = rownames(seu_demux[["HTO"]])[1:2],
  ncol = 3
)

colnames(seu_demux@meta.data)


# subset singlet only
seu_singlet <- subset(seu_demux, idents = c("tumor", "csf"))
head(seu_singlet@meta.data)

# subset csf and tumor samples into dif. variables (for tcr analysis)
seu_singlet_csf <- subset(seu_singlet, subset = MULTI_ID == "csf")
seu_singlet_tumor <- subset(seu_singlet, subset = MULTI_ID == "tumor")

#
# analysis of tumor and csf together!!!!
DefaultAssay(seu_singlet) <- "RNA"
seu_singlet <- NormalizeData(seu_singlet)
seu_singlet <- FindVariableFeatures(seu_singlet)
seu_singlet <- ScaleData(seu_singlet)

seu_singlet <- RunPCA(seu_singlet)
seu_singlet <- FindNeighbors(
  seu_singlet,
  reduction = "pca", dims = 1:50
)
seu_singlet <- FindClusters(
  seu_singlet,
  resolution = 1.1, verbose = FALSE
)
seu_singlet <- RunUMAP(
  seu_singlet,
  reduction = "pca", dims = 1:30
)

DimPlot(
  seu_singlet,
  group.by = "MULTI_ID", reduction = "umap"
)


library(viridis)
mycolor2 <- scale_color_viridis(option = "viridis")

FeaturePlot(
  seu_singlet, c("CD8B", "NKG7"),
  max.cutoff = 2
) & mycolor2

FeaturePlot(
  seu_singlet, c("CD4"),
  max.cutoff = 2
) & mycolor2

FeaturePlot(
  seu_singlet, c("CX3CR1"),
  max.cutoff = 2
) & mycolor2
