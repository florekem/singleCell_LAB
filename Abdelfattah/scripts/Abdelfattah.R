# Chromium 3` library

library(Seurat)
library(ggplot2)
library(Matrix)
library(tidyverse)
library(clustree)
library(scCustomize)
library(parallel)

setwd("/mnt/sda4/singleCell_LAB/Abdelfattah")

s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

seu_list <- list()
sample_meta_file <- read.csv("scripts/file_names.csv")
sample_names <- sample_meta_file$name

for (i in sample_names) {
  print(i)
  sample_path <- file.path("data", i)
  print(sample_path)
  sparse_matrix <- Seurat::Read10X(data.dir = sample_path)
  seu_list[[i]] <- CreateSeuratObject(
    counts = sparse_matrix,
    project = paste0(i),
    min.cells = 5,
    min.features = 300
  )
  if (grepl("rGBM", i)) {
    seu_list[[i]]$glioma_type <- "rGBM"
  }
  if (grepl("ndGBM", i)) {
    seu_list[[i]]$glioma_type <- "ndGBM"
  }
  if (grepl("LGG", i)) {
    seu_list[[i]]$glioma_type <- "LGG"
  }
}

# merging based on:
# https://github.com/samuel-marsh/scCustomize/blob/HEAD/R/Object_Utilities.R
# first, rename cells, bc it wont be possible using reduce
seu_list <- lapply(seq_along(seu_list), function(x) {
  seu_list[[x]] <- RenameCells(
    object = seu_list[[x]],
    add.cell.id = sample_meta_file$name[x]
  )
})
# merge list using reduce function
# Reduce a list to a single value by iteratively applying a binary function
# as I understand, reduce is doing next() on elements of list
seu_merged <- reduce(seu_list, function(x, y) {
  merge(
    x = x,
    y = y,
    merge.data = merge.data, # TRUE
    project = "merged",
  )
})
# looks the same as Merge_Seurat_List() from scCustomize package,
# but it skips initial checks implemeted there.

seu_merged
head(rownames(seu_merged@meta.data))
DefaultAssay(seu_merged)

seu_merged <- PercentageFeatureSet(
  seu_merged,
  pattern = "^MT-",
  col.name = "percent_mito"
)

seu_merged <- PercentageFeatureSet(
  seu_merged,
  pattern = "^RP[SL]",
  col.name = "percent.ribo"
)

seu_merged[["log10GenesPerUmi"]] <-
  log10(seu_merged$nFeature_RNA) / log10(seu_merged$nCount_RNA)

colnames(seu_merged@meta.data)

VlnPlot(
  seu_merged,
  features = c(
    "nFeature_RNA", "nCount_RNA",
    "percent_mito", "percent.ribo", "log10GenesPerUmi"
  )
)

seu_merged <- subset(
  seu_merged,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 6000 &
    percent_mito < 20 &
    percent.ribo < 50
  # log10umi?
)

# seu_merged <- JoinLayers(seu_merged)

# doubletfinder

seu_merged <- NormalizeData(seu_merged)

# seu_merged <- CellCycleScoring(
#   object = seu_merged,
#   s.features = s_genes,
#   g2m.features = g2m_genes,
#   set.ident = TRUE
# )

# seu_merged$CC.Difference <- seu_merged$S.Score - seu_merged$G2M.Score

seu_merged <- FindVariableFeatures(
  seu_merged,
  selection.method = "vst",
  n.Features = 2000
)

seu_merged <- ScaleData(seu_merged)

seu_merged <- RunPCA(seu_merged)

seu_merged <- IntegrateLayers(
  seu_merged,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony"
)

# saveRDS(seu_merged, file = "rds/harmony_integrated.rds")

seu_merged <- readRDS("rds/harmony_integrated.rds")

seu_merged_list <- list()
dims <- list(10, 20, 30)

# multicore do not work, too low ram, but why?
# single-core do not use as much memory ~22gb
seu_merged_list <- lapply(seq_along(dims), function(x) {
  seu_merged_list[[x]] <- Seurat::FindNeighbors(
    seu_merged,
    reduction = "pca", dims = 1:dims[[x]]
  )
})

seu_merged <- FindNeighbors(
  seu_merged,
  reduction = "pca", dims = 1:10
)

seu_merged_list <- lapply(seq_along(dims), function(x) {
  seu_merged_list[[x]] <- FindClusters(
    seu_merged_list[[x]],
    resolution = 0.3, verbose = FALSE
  )
})

seu_merged_list <- lapply(seq_along(dims), function(x) {
  RunUMAP(
    seu_merged_list[[x]],
    dims = 1:dims[[x]]
  )
})


DimPlot(
  seu_merged_list[[3]],
  label = TRUE,
  group.by = "RNA_snn_res.0.3"
)



FeaturePlot(
  seu_merged,
  features = "percent_mito",
  label = TRUE
)


seu_rgbm <- subset(seu_merged, subset = glioma_type == "rGBM")
seu_rgbm

seu_rgbm <- FindVariableFeatures(
  seu_rgbm,
  selection.method = "vst",
  n.Features = 2000
)

seu_rgbm <- ScaleData(seu_rgbm)

seu_rgbm <- RunPCA(seu_rgbm)

seu_rgbm <- RunUMAP(
  seu_rgbm,
  dims = 1:30
)
DimPlot(seu_rgbm)



# https://github.com/chris-mcginnis-ucsf/DoubletFinder
# Ensure that input data is cleared of low-quality cell clusters.
# There are a variety of ways to do this, but I usually use the
# following workflow:
# 1. Manually threshold raw gene expression matrices according to RNA nUMIs
# (especially important when dealing with super-loaded 10X data
# because of the way CellRanger threholds data -- See Lun et al., 2019, Genome Biology.
# 2. Pre-process data using standard workflow.
# Identify clusters with (A) low RNA UMIs,
# (B) High % mitochondrial reads, and/or
# (C) Uninformative marker genes.
# Remove clusters, pre-process again, and run DoubletFinder.
