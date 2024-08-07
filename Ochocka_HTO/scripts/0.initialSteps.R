# InitialSteps ------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(Matrix)
library(tidyverse)
library(clustree)

setwd("/mnt/sda4/singleCell_LAB/Ochocka_HTO")

# prepare list of Seurat objects.
seu_list <- list()
sample_names <- c("replicate1", "replicate2", "replicate3")

for (i in sample_names) {
  print(i)
  sample_path <- file.path("data", i)
  print(sample_path)
  sparse_matrix <- Seurat::Read10X(data.dir = sample_path)
  gex <- sparse_matrix$`Gene Expression`
  adt <- sparse_matrix$`Antibody Capture`[grepl(
    "TotalSeqB", rownames(sparse_matrix$`Antibody Capture`)
  ), ]

  hto <- sparse_matrix$`Antibody Capture`[grepl(
    "Hashtag", rownames(sparse_matrix$`Antibody Capture`)
  ), ]
  rownames(hto) <- c("F_ctrl", "M_ctrl", "F_D14", "M_D14", "F_D21", "M_D21")

  joint_bcs <- intersect(colnames(gex), colnames(hto))

  gex <- gex[, joint_bcs]
  adt <- adt[, joint_bcs]
  hto <- hto[, joint_bcs]
  seu_list[[i]] <- CreateSeuratObject(
    counts = gex,
    project = paste0(i),
    min.cells = 5,
    min.features = 300
  )
  hto <- hto[, rownames(seu_list[[i]]@meta.data)]
  seu_list[[i]][["HTO"]] <- CreateAssay5Object(counts = hto)
  adt <- adt[, rownames(seu_list[[i]]@meta.data)]
  seu_list[[i]][["ADT"]] <- CreateAssay5Object(counts = adt)
}

seu_merged <- merge(
  seu_list[[1]], c(seu_list[[2]], seu_list[[3]]),
  add.cell.ids = c("replicate1", "replicate2", "replicate3")
)

head(rownames(seu_merged@meta.data))
# "replicate1_ACCACAAAGAGATGCC-1"
seu_merged
# An object of class Seurat
# 18333 features across 36720 samples within 3 assays
# Active assay: HTO (6 features, 0 variable features)
# 6 layers present: counts.replicate1, counts.replicate2,
# counts.replicate3, data.replicate1, data.replicate2, data.replicate3
# 2 other assays present: RNA, ADT

DefaultAssay(seu_merged) <- "HTO"

seu_merged <- PercentageFeatureSet(
  seu_merged,
  pattern = "^mt-",
  col.name = "percent_mito"
)
seu_merged[["log10GenesPerUmi"]] <-
  log10(seu_merged$nFeature_RNA) / log10(seu_merged$nCount_RNA)

seu_joined <- JoinLayers(seu_merged)
seu_joined <- NormalizeData(
  seu_joined,
  assay = "HTO", normalization.method = "CLR"
)
seu_joined <- HTODemux(
  seu_joined,
  assay = "HTO",
  positive.quantile = 0.99
)
seu_joined

RidgePlot(
  seu_joined,
  assay = "HTO",
  features = rownames(seu_joined[["HTO"]])[1:6],
  ncol = 3
)

seu_joined$sex <- "test" # need to create filled columns first
seu_joined$condition <- "test"
seu_joined$condition[grepl("ctrl", seu_joined$hash.ID)] <- "ctrl"
seu_joined$condition[grepl("D14", seu_joined$hash.ID)] <- "D14"
seu_joined$condition[grepl("D21", seu_joined$hash.ID)] <- "D21"
seu_joined$sex[grepl("F", seu_joined$hash.ID)] <- "Female"
seu_joined$sex[grepl("M", seu_joined$hash.ID)] <- "Male"
colnames(seu_joined@meta.data)
head(seu_joined@meta.data)

HTOHeatmap(seu_joined, assay = "HTO", ncells = 5000)

# subset singlets only
Idents(seu_joined) <- "HTO_classification.global"

seu_joined_singlet <- subset(seu_joined, idents = "Singlet")

RidgePlot(
  seu_joined_singlet,
  assay = "HTO",
  features = rownames(seu_joined_singlet[["HTO"]])[1:6],
  ncol = 3
)

# Integrate RNA ---------------
# jako że joinowałem tylko assay hto, mogę bezkarnie integrować, bez dzielenia
DefaultAssay(seu_joined_singlet) <- "RNA"
seu_joined_singlet <- NormalizeData(seu_joined_singlet)
seu_joined_singlet <- FindVariableFeatures(seu_joined_singlet)
seu_joined_singlet <- ScaleData(seu_joined_singlet)

seu_joined_singlet <- RunPCA(seu_joined_singlet)
seu_joined_singlet <- FindNeighbors(
  seu_joined_singlet,
  reduction = "pca", dims = 1:50
)
seu_joined_singlet <- FindClusters(
  seu_joined_singlet,
  resolution = 1.1, verbose = FALSE
)
seu_joined_singlet <- RunUMAP(
  seu_joined_singlet,
  reduction = "pca", dims = 1:30
)

seu_joined_singlet <- IntegrateLayers(
  object = seu_joined_singlet, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  assay = "RNA",
  verbose = FALSE
)
seu_joined_singlet <- FindNeighbors(
  seu_joined_singlet,
  reduction = "integrated.cca", dims = 1:30
)
seu_joined_singlet <- FindClusters(
  seu_joined_singlet,
  resolution = 1.1, verbose = FALSE, cluster.name = "cca_clusters"
)
seu_joined_singlet <- RunUMAP(
  seu_joined_singlet,
  reduction = "integrated.cca", dims = 1:10, reduction.name = "umpa.cca"
)

DimPlot(
  seu_joined_singlet,
  group.by = "orig.ident", reduction = "umpa.cca"
)
DimPlot(
  seu_joined_singlet,
  group.by = "seurat_clusters", split.by = "condition", reduction = "umpa.cca"
)
DimPlot(
  seu_joined_singlet,
  group.by = "seurat_clusters", reduction = "umpa.cca"
)

# Visualize multiple modalities side-by-side ----------
DefaultAssay(seu_joined_singlet) <- "ADT"
seu_joined_singlet <- JoinLayers(seu_joined_singlet)
seu_joined_singlet <- NormalizeData(
  seu_joined_singlet,
  normalization.method = "CLR", margin = 2
)

rownames(seu_joined_singlet)
DefaultAssay(seu_joined_singlet) <- "ADT"
rownames(seu_joined_singlet)

p1 <- FeaturePlot(
  seu_joined_singlet, "CD49d-TotalSeqB",
  max.cutoff = 2
) & mycolor2

DefaultAssay(seu_joined_singlet) <- "RNA"
p2 <- FeaturePlot(seu_joined_singlet, "Itga4", max.cutof = 2, cols = mycolor)

p1 | p2


library(viridis)
mycolor2 <- scale_color_viridis(option = "viridis")

mycolor <- colorRampPalette(rev(brewer.pal(n = 5, name = "RdYlBu")))(100)
mycolor <- rev(brewer.pal(n = 10, name = "Spectral"))
