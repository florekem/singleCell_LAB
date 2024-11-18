library(Seurat)
library(ggplot2)
library(Matrix)
library(tidyverse)
library(clustree)

setwd("/mnt/sda4/singleCell_LAB/HitGlio_sample1-test")

# prepare list of Seurat objects.
# seu_list <- list()
# sample_names <- c("sample1")

sample_path <- file.path("data", "sample1")
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
tail(hto)
adt

head(gex)
dim(gex)
dim(adt)
dim(hto)

joint_bcs <- intersect(colnames(gex), colnames(hto))
joint_bcs
gex <- gex[, joint_bcs]
ncol(gex)
adt <- adt[, joint_bcs]
ncol(adt)
hto <- hto[, joint_bcs]
ncol(hto)

seu <- CreateSeuratObject(
  counts = gex,
  project = "sample1",
  min.cells = 10,
  min.features = 30
)
hto <- hto[, rownames(seu@meta.data)]
seu[["HTO"]] <- CreateAssay5Object(counts = hto)
adt <- adt[, rownames(seu@meta.data)]
seu[["ADT"]] <- CreateAssay5Object(counts = adt)

# seu_merged <- merge(
# seu_list[[1]], c(seu_list[[2]], seu_list[[3]]),
# add.cell.ids = c("replicate1", "replicate2", "replicate3")
# )

head(rownames(seu@meta.data))
DefaultAssay(seu) <- "HTO"
seu
# seu <- PercentageFeatureSet(
#   seu,
#   pattern = "^mt-",
#   col.name = "percent_mito"
# )
seu[["log10GenesPerUmi"]] <-
  log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)

# seu_joined <- JoinLayers(seu)

seu <- NormalizeData(
  seu,
  assay = "HTO", normalization.method = "CLR"
)
seu <- HTODemux(
  seu,
  assay = "HTO",
  positive.quantile = 0.99
)
seu@meta.data

RidgePlot(
  seu,
  assay = "HTO",
  features = rownames(seu[["HTO"]])[1:2],
  ncol = 3
)
table(seu$hash.ID)
HTOHeatmap(seu, assay = "HTO", ncells = 5000)

# subset singlets only
Idents(seu) <- "HTO_classification.global"
table(seu$HTO_classification.global)

VlnPlot(seu, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

seu_singlet <- subset(seu, idents = "Singlet")

RidgePlot(
  seu_singlet,
  assay = "HTO",
  features = rownames(seu_singlet[["HTO"]])[1:2],
  ncol = 1
)

HTOHeatmap(seu_singlet, assay = "HTO", ncells = 5000)

seu_singlet
seu
