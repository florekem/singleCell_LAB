# according to: https://satijalab.org/seurat/articles/hashing_vignette

# httpgd::hgd_url() #nolint

# 0. Libraries --------------------------------------------------
library(Seurat)
library(ggplot2)
library(Matrix)
library(tidyverse)
library(clustree)
library(reticulate) # load packages required by leiden algo
reticulate::use_condaenv("/mnt/sda4/conda/envs/R")
library(viridis) # colors

options(future.globals.maxSize = 2.0 * 1e9) # 2GB
library(sceasy) # convert to anndata, scvi integration
library(dittoSeq)
# 1. Paths ------------------------------------------------------
setwd("/mnt/sda4/singleCell_LAB/HitGlio_sample1-test")

sample_path <- file.path("/mnt/sdb1/runs/sample1_multilane_spec_index_NNN/outs/per_sample_outs/sample1/count/sample_filtered_feature_bc_matrix")
print(sample_path)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

source("scripts/functions.R")

# Podejscie 1
seu <- initial_steps(sample_path)
seu_demux <- initial_demux(seu)
Idents(seu_demux) <- "MULTI_ID"

VlnPlot(seu_demux, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# head(seu_demux@meta.data)
# seu_singlet <- seu_subset(
#   seu_demux,
#   idents = c("tumor", "csf")
# )
seu_singlet <- subset(seu_demux, idents = c("tumor", "csf"))

head(seu_singlet@meta.data)

### subset csf and tumor singlets (for tcr analysis) -------
seu_singlet_csf <- seu_subset(
  seu_singlet,
  subset = MULTI_ID == "csf"
)
seu_singlet_tumor <- seu_subset(
  seu_singlet,
  subset = MULTI_ID == "tumor"
)
# /now this goes to the script where tcr data in analyzed/

seu_singlet <- add_quality_features(seu_singlet)
colnames(seu_singlet@meta.data)

# może warto by wyczyścić scale.data, żeby nie wyświetlało warning?
# temp. bo nie bede potrzebowal tego pozniej
seu_singlet_q_filt_temp <- seu_normalize_var_scale_pca(
  seu_singlet,
  normalize = "log",
  cluster.name = "cluster.log",
  reduction = "pca",
  reduction.name = "umap.log",
  dims = 1:30,
  k.param = 20,
  algorithm = 4,
  resolution = 0.5,
  group.singletons = TRUE
)
seu_singlet_q_filt_temp

ElbowPlot(seu, ndims = 40)

## Visualize clusters and features before quality filtering ------------------
q_features <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "percent_mito",
  "percent_ribo",
  "percent_hemoglob"
)

p2 <- seu_DimPlot(
  seu_singlet_q_filt_temp,
  group.by = "MULTI_ID",
  reduction = "umap"
)

p3 <- seu_FeaturePlot(
  seu_singlet_q_filt_temp,
  q_features,
  reduction = "umap"
)

p4 <- violin_plot(
  seu_singlet_q_filt_temp,
  features = q_features
)

## 5.5. Subset cells based on quality assesment -------------
# subset juz oryginalnego singlet
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
seu_singlet_q_filt_temp <- seu_normalize_var_scale_pca(
  seu_singlet,
  normalize = "log",
  cluster.name = "cluster.log",
  reduction = "pca",
  reduction.name = "umap.log",
  dims = 1:30,
  k.param = 20,
  algorithm = 4,
  resolution = 0.5,
  group.singletons = TRUE
)

seu_singlet <- NormalizeData(
  seu_singlet,
  assay = "ADT",
  normalization.method = "CLR"
)


## 6.1 Save RDS after quality filtering -------------------------------
saveRDS(seu_singlet, "RDS/seu_singlet_after-filtering.rds")

## 6.2. Vizualize clusters and quality features after quality filtiring -----

# 7. Use SCT transform (replace 6.) ------------
seu_singlet <- seu_normalize_var_scale_pca(
  seu_singlet,
  normalize = "sct",
  cluster.name = "cluster.sct",
  pca.reduction.name = "pca.sct",
  umap.reduction.name = "umap.sct",
  dims = 1:30,
  k.param = 30,
  algorithm = 4,
  resolution = 0.5,
  group.singletons = FALSE
)
seu_singlet
?RunPCA
seu_singlet <- NormalizeData(
  seu_singlet,
  assay = "ADT",
  normalization.method = "CLR"
)

## 7.3 Save RDS of clustered singlets
saveRDS(seu_singlet, "RDS/seu_singlet-clustered-27nov2024.rds")
seu_singlet

# 8. Visualize GEX and ADT features ---------
## 8.1 GEX ----------------
DefaultAssay(seu_singlet) <- "SCT"
inferno <- viridis::scale_color_viridis(option = "inferno")

vars <- c("CD8A", "CD4")
seu_FeaturePlot(
  seu_singlet,
  assay = "SCT",
  features = vars,
  reduction = "umap.sct",
  color = inferno,
  max.cutoff = NA,
  show = TRUE,
  save_path = FALSE
)
seu_FeaturePlot(seu_singlet, "SCT", vars)

seu_dittoDotPlot(
  seu_singlet,
  vars = vars,
  group.by = "seurat_clusters"
)

## ADT -------------
features <- Features(seu_singlet)
length(features)
f1 <- features[1:10]
f2 <- features[11:21]
f3 <- features[22:32]
f4 <- features[33:42]

seu_FeaturePlot(
  seu_singlet, "ADT",
  reduction = "umap.sct",
  color = inferno,
  "TotalSeqC-CD4"
)


# 9. totalVI model training ------------------
# /not sure if it will produce any good output, but lets try

## 9.1 convert to anndata (mudata) -----------
DefaultAssay(seu_singlet) <- "SCT"

## for the purpose of scVI extract 3k variable genes (stored in scale.data)
## and add it to the seurat object (TRUE/FALSE column)
scaled_matrix <- GetAssayData(seu_singlet, layer = "scale.data", assay = "SCT")
test <- rownames(seu_singlet) %in% rownames(scaled_matrix)
seu_singlet <- AddMetaData(seu_singlet,
  metadata = test,
  col.name = "sct.var.genes"
)
seu_singlet@meta.data$sct.var.genes

# seu_singlet <- readRDS("RDS/seu_singlet-clustered-27nov2024.rds")

# seu_singlet[["RNA3"]] <- as(object = seu_singlet[["RNA"]], Class = "Assay")
# DefaultAssay(seu_singlet) <- "RNA3"
# seu_singlet[["RNA"]] <- NULL
# seu_singlet <- RenameAssays(object = seu_singlet, RNA3 = "RNA")

DefaultAssay(seu_singlet) <- "ADT"
seu_singlet
adata <- convertFormat(
  seu_singlet,
  from = "seurat",
  to = "anndata",
  main_layer = "counts",
  # assay  = "SCT",
  assay = "ADT",
  drop_single_values = FALSE,
  # outFile = "RDS/test-adt.h5ad"
  outFile = "RDS/scvi-adt_12-02-24.h5ad"
)
adata

## go to python, and bring back latent.csv
latent <- read.csv("RDS/latent-sct-totalVI_12-02-24.csv")
latent_mtx <- as.matrix(latent)
rownames(latent_mtx) <- colnames(seu_singlet)
latent_mtx <- latent_mtx[, -1] # remove 1st col
dim(latent_mtx)
colnames(latent_mtx) <- paste0("totalvi_", 1:20)
colnames(latent_mtx)

DefaultAssay(seu_singlet) <- "SCT"
seu_singlet
colnames(seu_singlet@meta.data)

seu_singlet[["totalvi"]] <- CreateDimReducObject(
  embeddings = latent_mtx,
  key = "totalvi_",
  assay = DefaultAssay(seu_singlet)
)
seu_singlet
seu_singlet <- seu_normalize_var_scale_pca(
  seu_singlet,
  normalize = FALSE,
  cluste.name = "cluster.totalvi",
  reduction = "totalvi",
  reduction.name = "umap.totalvi",
  dims = 1:20,
  k.param = 30,
  algorithm = 4,
  resolution = 0.8,
  group.singletons = FALSE
)

# VIZ. OF scVI data
DimPlot(
  seu_singlet,
  reduction = "umap",
  pt.size = 0.5,
  group.by = c("MULTI_ID", "seurat_clusters"),
  label = TRUE
)

DimPlot(
  seu_singlet,
  reduction = "umap",
  pt.size = 0.5,
  group.by = "signleR_anot",
  label = TRUE
)

DefaultAssay(seu_singlet) <- "SCT"

FeaturePlot(
  seu_singlet,
  # tcels,
  macrophages,
  reduction = "umap",
  label = TRUE
)

DefaultAssay(seu_singlet) <- "ADT"

seu_singlet <- NormalizeData(
  seu_singlet,
  assay = "ADT",
  normalization.method = "CLR"
)

FeaturePlot(
  seu_singlet,
  # c("TotalSeqC-CD8", "TotalSeqC-CD4"),
  c("TotalSeqC-CD8", "TotalSeqC-CD14"),
  # max.cutoff = 4,
  label = TRUE
) & mycolor2

# 10. vizualize clonality ----------------------------------------
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


# 11. Name clusters -------------------------------------
## copy seurat_clusters column
seu_singlet[["my_id"]] <- seu_singlet[["seurat_clusters"]]
Idents(seu_singlet) <- "my_id"
seu_singlet$my_id <- plyr::mapvalues(
  Idents(seu_singlet),
  from = c(1:23),
  to = c(
    "CD8",
    "CD4",
    "3",
    "CD11c myelo",
    "5",
    "NK",
    "CD4",
    "act. microglia",
    "monocytes/neutrophils",
    "CD8",
    "CD4",
    "12",
    "CD4",
    "CD4",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "23"
  )
)

DimPlot(
  seu_singlet,
  reduction = "umap",
  pt.size = 0.5,
  group.by = c("MULTI_ID", "my_id"),
  label = TRUE
)
colnames(seu_singlet@meta.data)

# Cell cycle scoring ------------------------------
# should be performed much earlier.
# as i dont think i should regress out phase genes,
# i left it here for the future (or mayby i'll change my mind)
seu_singlet <- CellCycleScoring(
  seu_singlet,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)

FeaturePlot(
  seu_singlet,
  c("S.Score", "G2M.Score"),
  reduction = "umap",
  label = TRUE
)

test <- RunPCA(seu_singlet, features = c(s.genes, g2m.genes))
DimPlot(test, reduction = "pca")

DimPlot(test,
  reduction = "pca",
  group.by = "Phase",
  split.by = "Phase"
)

FeaturePlot(
  seu_singlet, "percent_ribo",
  reduction = "umap",
  label = TRUE
)

Idents(seu_singlet) <- "my_id"
FeaturePlot(
  seu_singlet,
  # c("CX3CR1", "P2RY12"),
  # c("IL1B", "TNF", "CCL4"),
  # c("GFAP"),
  # c("FOLR2", "LYVE1", "MRC1", "CD206"),
  # "CD8B",
  tcels,
  # "PDCD1",
  # "CADM1",
  # c("TREM1", "TREM2"),
  # c("p53", "ATM", "CHEK2"),
  # c("HIF1A"),
  # c("BAX", "CASP3", "FAS"),
  # c("MICA", "MICB"),
  # c("TREM2", "SPP1"),
  # c("CD34", "CD38", "CD123", "CD117", "CD13", "CD33", "HLA‐DR"),
  reduction = "umap",
  label = TRUE
)

cluster_microglia.markers <- FindMarkers(
  seu_singlet,
  ident.1 = "microglia",
  min.pct = 0.25
)
rownames(cluster_microglia.markers)[1:30]

cluster_3.markers <- FindMarkers(
  seu_singlet,
  ident.1 = "3",
  min.pct = 0.25
)
rownames(cluster_3.markers)[1:30]

saveRDS(seu_singlet, "RDS/seu_singlet-clustered-totalvi-28nov2024.rds")
