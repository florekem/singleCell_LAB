# https://www.10xgenomics.com/datasets/20-k-bone-marrow-mononuclear-cells-bmmn-cs-5-ht-v-2-0-2-high-6-1-0 #nolint

# variable names with "_", variable_name
# column names with "." column.name

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
colnames(seu@meta.data)
# --- filter seuratObject ----------------------------------------------------
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

seu <- PercentageFeatureSet(
  seu,
  pattern = "^MT-",
  col.name = "percent.mito"
)

seu <- PercentageFeatureSet(
  seu,
  pattern = "^RP[SL]",
  col.name = "percent.ribo"
)

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo")
)

seu <- subset( # not sure about mito percentage
  seu,
  subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mito < 10
) # ew. jeszcze ncount > 20k
seu
# An object of class Seurat
# 36601 features across 16676 samples within 1 assay
# Active assay: RNA (36601 features, 0 variable features)
#  1 layer present: counts
VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo")
)

seu <- NormalizeData(seu)

seu <- CellCycleScoring(
  object = seu,
  s.features = s_genes,
  g2m.features = g2m_genes,
  set.ident = TRUE
)

seu$CC.Difference <- seu$S.Score - seu$G2M.Score

options(future.globals.maxSize = 891289600)
seu <- SCTransform(
  seu,
  verbose = FALSE,
  vars.to.regress = c(
    "nCount_RNA",
    "percent.mito",
    "CC.Difference" # not sure is ok to do that,
    # but it is usefull for Merging cycling
    # and non-cycling cells of the same type in one cluster
    # https://sib-swiss.github.io/single-cell-training/day2/day2-4_cell_annotation.html # nolint
  )
)
seu
# An object of class Seurat
# 55862 features across 16676 samples within 2 assays
# Active assay: SCT (19261 features, 3000 variable features)
#  3 layers present: counts, data, scale.data
#  1 other assay present: RNA
seu <- RunPCA(seu)

Seurat::ElbowPlot(seu, ndims = 40)
Seurat::DimHeatmap(seu, dims = 1:10, cells = 500, balanced = TRUE)

seu <- FindNeighbors(
  seu,
  dims = 1:30, reduction = "pca"
)

seu <- FindClusters(
  seu,
  algorithm = 1,
  # algorithm = 4, # problem with installing laiden for now
  resolution = seq(0.1, 1.8, by = 0.1)
)

clustree::clustree(
  seu@meta.data[, grep("SCT_snn_res", colnames(seu@meta.data))],
  prefix = "SCT_snn_res."
)

seu <- RunUMAP(
  object = seu,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap"
)

DimPlot(seu,
  reduction = "umap",
  group.by = "SCT_snn_res.0.3",
  # split.by = c("condition"),
  label = TRUE
)

tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
FeaturePlot(seu, tcell_genes, ncol = 2)

monocyte_genes <- c("CD14", "CST3", "CD68", "CTSS")
FeaturePlot(seu, monocyte_genes, ncol = 2)

DimPlot(seu, group.by = "Phase") # cclycle viz.

ref <- celldex::NovershternHematopoieticData()
class(ref)
table(ref$label.main)

seu_SingleR <- SingleR::SingleR(
  test = Seurat::GetAssayData(seu),
  ref = ref,
  labels = ref$label.main
)
dev.off()
SingleR::plotScoreHeatmap(seu_SingleR)

singleR_labels <- seu_SingleR$labels
t <- table(singleR_labels)
other <- names(t)[t < 10]
singleR_labels[singleR_labels %in% other] <- "none"
seu$SingleR_annot <- singleR_labels

dittoSeq::dittoDimPlot(seu, "SingleR_annot", size = 0.7)
dittoSeq::dittoBarPlot(seu, var = "SingleR_annot", group.by = "orig.ident")
dittoSeq::dittoBarPlot(seu, var = "SingleR_annot", group.by = "SCT_snn_res.0.3")



# --- scRep prepare T data -----------------------------------------------------
tcr_data_list <- read.csv("singleCell_LAB/scREP_projectTIL_tut/data/20k_BMMNCs/tcr/vdj_t_filtered_contig_annotations.csv") # nolint

combined_tcr <- combineTCR(
  tcr_data_list,
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

seu_tcr <- combineExpression(
  combined_tcr,
  seu,
  cloneCall = "gene",
  proportion = TRUE
)

colorblind_vector <- hcl.colors(n = 7, palette = "inferno", fixup = TRUE)
DimPlot(
  seu_tcr,
  group.by = "cloneSize"
) + scale_color_manual(values = rev(colorblind_vector[c(1, 3, 4, 5, 7)]))

# CD4
ref_cd4 <- load.reference.map(
  "singleCell_LAB/scREP_projectTIL_tut/data/refs/pTILs_hsa_cd4t.rds"
)

seu_tcr <- Run.ProjecTILs( # these variable names doesnt make any sense, fix it
  seu_tcr,
  ref = ref_cd4,
  ncores = 1
)

p1 <- plot.projection(ref_cd4) + theme(aspect.ratio = 1)
p2 <- plot.projection(
  ref_cd4, seu_tcr,
  linesize = 0.3, pointsize = 0.5
) + theme(aspect.ratio = 1)
p1 | p2

# Look at distribution of T cells in terms of cell states.
plot.statepred.composition(ref_cd4, seu_tcr, metric = "Percent")

# We can check the gene expression profile of the cells assigned
# to each state (yellow), and compare them to those of the
# reference states (black).
plot.states.radar(
  ref = ref_cd4, seu_tcr, min.cells = 30
)

# CD8
ref_cd8 <- load.reference.map(
  "singleCell_LAB/scREP_projectTIL_tut/data/refs/pTILs_hsa_cd8t.rds"
)

seu_tcr <- Run.ProjecTILs( # these variable names doesnt make any sense, fix it
  seu_tcr,
  ref = ref_cd8,
  ncores = 1
)

p1 <- plot.projection(ref_cd8) + theme(aspect.ratio = 1)
p2 <- plot.projection(
  ref_cd8, seu_tcr,
  linesize = 0.3, pointsize = 0.5
) + theme(aspect.ratio = 1)
p1 | p2

plot.statepred.composition(ref_cd8, seu_tcr, metric = "Percent")

plot.states.radar(
  ref = ref_cd8, seu_tcr, min.cells = 30
)

# --- scRep prepare B data -----------------------------------------------------
