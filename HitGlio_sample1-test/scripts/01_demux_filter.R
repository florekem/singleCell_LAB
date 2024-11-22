# according to: https://satijalab.org/seurat/articles/hashing_vignette

# httpgd::hgd_url() #nolint

# 0. Libraries --------------------------------------------------
library(Seurat)
library(ggplot2)
library(Matrix)
library(tidyverse)
library(clustree)
library(reticulate) # load packages required by leiden algo
library(leiden) # not sure if this one is req, check that next
use_python("/usr/bin/python3") # this is crucial for leiden working
library(viridis) # colors
mycolor2 <- scale_color_viridis(option = "viridis")
library(dittoSeq)
# 1. Paths ------------------------------------------------------
setwd("/mnt/sda4/singleCell_LAB/HitGlio_sample1-test")

sample_path <- file.path("/mnt/sdb1/runs/sample1_multilane_spec_index_NNN/outs/per_sample_outs/sample1/count/sample_filtered_feature_bc_matrix")
print(sample_path)

# 2. Create sparse matrix ---------------------------------------
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

# 3. Create seurat object ----------------------------------------
seu <- CreateSeuratObject(
  counts = Matrix::Matrix(as.matrix(gex), sparse = TRUE)
)

# 4. Demux and subset singlets -----------------------------------
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "mean.var.plot")
seu <- ScaleData(seu, features = VariableFeatures(seu))

## 4.1. Add HTO data as a new assay independent from RNA ----
seu[["HTO"]] <- CreateAssayObject(counts = hto)
seu[["ADT"]] <- CreateAssayObject(counts = adt)

## 4.2. Normalize HTO data, with centered log-ratio (CLR) transformation ----
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

## 4.3 subset singlet only ----
seu_singlet <- subset(seu_demux, idents = c("tumor", "csf"))
head(seu_singlet@meta.data)

### 4.3.1 separ. subset csf and tumor singlets (for tcr analysis) -------
seu_singlet_csf <- subset(seu_singlet, subset = MULTI_ID == "csf")
seu_singlet_tumor <- subset(seu_singlet, subset = MULTI_ID == "tumor")
# /now this goes to the script where tcr data in analyzed/

# 5. Quality filter cells ----------------------------------------
## 5.1. Add quality features -------
### 5.1.1. Ribosomal --------------
seu_singlet <- PercentageFeatureSet(
  seu_singlet,
  pattern = "^RP[SL]",
  col.name = "percent_ribo"
)
### 5.1.2. Mitochondrial --------
seu_singlet <- PercentageFeatureSet(
  seu_singlet,
  pattern = "^MT-",
  col.name = "percent_mito"
)
### 5.1.3. Hemoglobin -------
seu_singlet <- PercentageFeatureSet(
  seu_singlet,
  pattern = "^HB[^(P)]",
  col.name = "percent_hemoglob"
)
### 5.1.4. log10 genes per umi ------------------
seu_singlet[["log10GenesPerUmi"]] <-
  log10(seu_singlet$nFeature_RNA) / log10(seu_singlet$nCount_RNA)

## 5.2. Visualize quality features before filtering ---------------
Seurat::VlnPlot(
  seu_singlet,
  features = c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent_mito",
    "percent_ribo",
    "percent_hemoglob"
  )
)
## 5.3. Check how clusters looks like before quality filtering ----------
DefaultAssay(seu_singlet) <- "RNA"
seu_singlet <- NormalizeData(seu_singlet)
seu_singlet <- FindVariableFeatures(seu_singlet)
seu_singlet <- ScaleData(seu_singlet)
seu_singlet <- RunPCA(seu_singlet)
ElbowPlot(seu_singlet, ndims = 40)

seu_singlet <- FindNeighbors(
  seu_singlet,
  reduction = "pca", dims = 1:30
)
seu_singlet <- FindClusters(
  seu_singlet,
  algorithm = 4, # leiden
  resolution = 0.5, verbose = FALSE
)
seu_singlet <- RunUMAP(
  seu_singlet,
  reduction = "pca", dims = 1:30
)

### 5.4. Visualize cluster before quality filtering ------------------
DimPlot(
  seu_singlet,
  group.by = "MULTI_ID", reduction = "umap"
)
FeaturePlot(
  seu_singlet, "percent_mito",
  reduction = "umap",
  label = TRUE
)
FeaturePlot(
  seu_singlet, "log10GenesPerUmi",
  reduction = "umap",
  label = TRUE
)
FeaturePlot(
  seu_singlet, "percent_ribo",
  reduction = "umap",
  label = TRUE
)
FeaturePlot(
  seu_singlet, "percent_hemoglob",
  reduction = "umap",
  label = TRUE
)
Seurat::VlnPlot(
  seu_singlet,
  features = c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent_mito",
    "percent_ribo",
    "percent_hemoglob"
  )
)
## 5.5. Subset cells based on quality assesment -------------
seu_singlet <- subset(
  seu_singlet,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 7000 &
    nCount_RNA < 50000 &
    nCount_RNA > 30 &
    percent_mito < 8
)
seu_singlet

# 6. Recluster after quality subsetting --------------------------------
seu_singlet <- NormalizeData(seu_singlet)
seu_singlet <- FindVariableFeatures(seu_singlet)
seu_singlet <- ScaleData(seu_singlet)
seu_singlet <- RunPCA(seu_singlet)
ElbowPlot(seu_singlet, ndims = 40)

seu_singlet <- FindNeighbors(
  seu_singlet,
  reduction = "pca", dims = 1:30
)
seu_singlet <- FindClusters(
  seu_singlet,
  resolution = 0.5, verbose = FALSE,
  algorithm = 4
)
seu_singlet <- RunUMAP(
  seu_singlet,
  reduction = "pca", dims = 1:30
)

## 6.1. Vizualize clusters and quality features after quality filtiring -----
DimPlot(
  seu_singlet,
  group.by = "MULTI_ID", reduction = "umap",
  label = TRUE
)
DimPlot(
  seu_singlet,
  group.by = "seurat_clusters", reduction = "umap",
  label = TRUE
)
FeaturePlot(
  seu_singlet, "percent_mito",
  reduction = "umap",
  label = TRUE
)
FeaturePlot( # w/o explicit removal of low umis, this parameter improoveed
  seu_singlet, "log10GenesPerUmi",
  reduction = "umap",
  label = TRUE
)
FeaturePlot(
  seu_singlet, "percent_ribo",
  reduction = "umap",
  label = TRUE
)
FeaturePlot(
  seu_singlet, "percent_hemoglob",
  reduction = "umap",
  label = TRUE
)
Seurat::VlnPlot(
  seu_singlet,
  features = c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent_mito",
    "percent_ribo",
    "percent_hemoglob"
  )
)

# 7. Use SCT transform (replace 6.) ------------
options(future.globals.maxSize = 2.0 * 1e9) # 2GB
?SCTransform
seu_singlet <- SCTransform(
  seu_singlet,
  verbose = FALSE
)

DefaultAssay(seu_singlet) <- "SCT"

seu_singlet <- RunPCA(seu_singlet)
ElbowPlot(seu_singlet, ndims = 40)

seu_singlet <- FindNeighbors(
  seu_singlet,
  reduction = "pca", dims = 1:50 # increased bc of SCT
)
seu_singlet <- FindClusters(
  seu_singlet,
  resolution = 0.2, verbose = FALSE,
  algorithm = 4 # leiden
)
seu_singlet <- RunUMAP(
  seu_singlet,
  reduction = "pca", dims = 1:50
)

## 7.1. Vizualize clusters and quality features after quality filtiring ------
DimPlot(
  seu_singlet,
  group.by = "MULTI_ID", reduction = "umap",
  label = TRUE
)
DimPlot(
  seu_singlet,
  group.by = "seurat_clusters", reduction = "umap",
  label = TRUE
)
FeaturePlot(
  seu_singlet, "percent_mito",
  reduction = "umap",
  label = TRUE
)
FeaturePlot( # w/o explicit removal of low umis, this parameter improoveed
  seu_singlet, "log10GenesPerUmi",
  reduction = "umap",
  label = TRUE
)
FeaturePlot(
  seu_singlet, "percent_ribo",
  reduction = "umap",
  label = TRUE
)
FeaturePlot(
  seu_singlet, "percent_hemoglob",
  reduction = "umap",
  label = TRUE
)
Seurat::VlnPlot(
  seu_singlet,
  features = c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent_mito",
    "percent_ribo",
    "percent_hemoglob"
  )
)

# 8. Visualize GEX and ADT features ---------
## 8.1 GEX ----------------
DefaultAssay(seu_singlet) <- "SCT"

active_gex <- f5
cutoff <- NA

FeaturePlot(
  seu_singlet,
  active_gex,
  max.cutoff = cutoff
) & mycolor2
?FeaturePlot

ditto_group_by <- "MULTI_ID"
ditto_group_by <- "seurat_clusters"
dittoSeq::dittoDotPlot(
  seu_singlet,
  vars = active_gex,
  group.by = ditto_group_by
)

## 8.2. ADT -------------
DefaultAssay(seu_singlet) <- "ADT"
DefaultAssay(seu_singlet)

seu_singlet <- NormalizeData(
  seu_singlet,
  assay = "ADT",
  normalization.method = "CLR"
)


features <- Features(seu_singlet)
length(features)
f1 <- features[1:10]
f2 <- features[11:21]
f3 <- features[22:32]
f4 <- features[33:42]

p1 <- FeaturePlot(
  seu_singlet, f1,
  # max.cutoff = 2,
  label = TRUE
) & mycolor2
p1
ggsave("plots/p4.png", p1)

saveRDS()

FeaturePlot(
  seu_singlet, "TotalSeqC-CD25",
  # max.cutoff = 4,
  label = TRUE
) & mycolor2





# vizualize clonality
library(scRepertoire)

tcr_sample <- read_csv("/mnt/sdb1/runs/sample1_multilane_spec_index_NNN/outs/per_sample_outs/sample1/vdj_t/filtered_contig_annotations.csv")

dim(tcr_sample)
# 16400 31

tcr_sample_list <- list(tcr_sample)

combined_tcr <- combineTCR(
  tcr_sample_list,
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE
)

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

seu_singlet <- combineExpression(
  combined_tcr,
  seu_singlet,
  cloneCall = "gene",
  proportion = TRUE
)

colorblind_vector <- hcl.colors(n = 7, palette = "inferno", fixup = TRUE)
DimPlot(
  seu_singlet,
  group.by = "cloneSize"
) + scale_color_manual(values = rev(colorblind_vector[c(1, 3, 4, 5, 7)]))
