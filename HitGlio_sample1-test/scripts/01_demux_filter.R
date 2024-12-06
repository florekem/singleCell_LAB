# according to: https://satijalab.org/seurat/articles/hashing_vignette

# httpgd::hgd_url() #nolint

source("scripts/functions.R")

# 0. Libraries --------------------------------------------------
library(Seurat)
library(ggplot2)
library(Matrix)
library(tidyverse)
library(clustree)
library(reticulate) # load packages required by leiden algo
reticulate::use_condaenv("/mnt/sda4/conda/envs/R")
library(viridis) # colors
library(DoubletFinder)
options(future.globals.maxSize = 2.0 * 1e9) # 2GB
library(sceasy) # convert to anndata, scvi integration
library(dittoSeq)

# 1. Paths ------------------------------------------------------
setwd("/mnt/sda4/singleCell_LAB/HitGlio_sample1-test")

sample_path <- file.path("/mnt/sdb1/runs/sample1_multilane_spec_index_NNN/outs/per_sample_outs/sample1/count/sample_filtered_feature_bc_matrix") # nolint
print(sample_path)

s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

# 2. Load Seurat Object, Demux ---------------------
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

# 3. Subset csf and tumor singlets (for tcr analysis) -------
## / this goes to the script where tcr data in analyzed/
seu_singlet_csf <- seu_subset(
  seu_singlet,
  subset = MULTI_ID == "csf"
)
seu_singlet_tumor <- seu_subset(
  seu_singlet,
  subset = MULTI_ID == "tumor"
)

# 4. Add quality features --------------------------
seu_singlet <- add_quality_features(seu_singlet)

# 5. Add Cell cycle scoring ------------------------------
seu_singlet <- CellCycleScoring(
  seu_singlet,
  s.features = s_genes,
  g2m.features = g2m_genes,
  set.ident = TRUE
)
colnames(seu_singlet@meta.data)

# 6. Normalize befeore quality filtering (check how it looks like)
# / może od razu SCT?
# / _temp bo nie bede potrzebowal tego pozniej
seu_singlet_q_filt_temp <- seu_normalize_var_scale_pca(
  seu_singlet,
  normalize = "log",
  cluster.name = "cluster.log",
  pca.reduction.name = "pca",
  umap.reduction.name = "umap.log",
  dims = 1:30,
  k.param = 20,
  algorithm = 4,
  resolution = 0.5,
  group.singletons = TRUE
)
seu_singlet_q_filt_temp

ElbowPlot(seu, ndims = 40)

# 7. Vizualize clusters and features before quality filtering ------------------
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

# 8. Subset cells based on quality assesment -------------
# subset juz oryginalnego singlet
seu_singlet <- subset(
  seu_singlet,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 7000 &
    nCount_RNA < 50000 &
    # nCount_RNA > 30 &
    nCount_RNA > 500 &
    percent_mito < 8
)
seu_singlet

# 9. Normalize using SCT after q filtering for DOUBLET REMOVAL ------------
seu_singlet <- seu_normalize_var_scale_pca(
  seu_singlet,
  normalize = "sct",
  cluster.name = "cluster.sct",
  pca.reduction.name = "pca.sct",
  umap.reduction.name = "umap.sct",
  dims = 1:10, # change to 10 in case of removing doublets
  k.param = 20,
  algorithm = 4,
  resolution = 0.8,
  group.singletons = TRUE
)
seu_singlet

seu_singlet <- NormalizeData(
  seu_singlet,
  assay = "ADT",
  normalization.method = "CLR"
)

# 7. Remove doublets based on SCT initial clustering --------------------------
homotypic_prop <- modelHomotypic(seu_singlet$cluster.sct)
nExp_poi <- round(0.075 * nrow(seu_singlet)) # estimate ~7.5% doublets
nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))

seu_singlet_temp <- DoubletFinder::doubletFinder(
  seu_singlet,
  PCs = 1:10,
  pN = 0.25,
  pK = 0.09,
  nExp = nExp_poi_adj,
  reuse.pANN = FALSE,
  sct = TRUE
)
seu_singlet_temp
colnames(seu_singlet_temp@meta.data)
head(seu_singlet_temp@meta.data)

seu_singlet <- subset(
  seu_singlet_temp,
  subset = DF.classifications_0.25_0.09_1841 == "Singlet"
)
seu_singlet
## / zostało 13409 samples z 15250

seu_singlet <- seu_normalize_var_scale_pca(
  seu_singlet,
  normalize = "sct",
  cluster.name = "cluster.sct",
  pca.reduction.name = "pca.sct",
  umap.reduction.name = "umap.sct",
  dims = 1:20,
  k.param = 40,
  algorithm = 4,
  resolution = 0.1,
  group.singletons = TRUE
)
seu_singlet <- NormalizeData(
  seu_singlet,
  assay = "ADT",
  normalization.method = "CLR"
)
seu_singlet

saveRDS(seu_singlet, "RDS/seu_singlet_sct-doublet-500_05dec24.rds")

# 8. totalVI model training of doublet-removed seu_singlet ------------------
DefaultAssay(seu_singlet) <- "SCT"

## /for the purpose of scVI extract 3k variable genes (stored in scale.data)
## /and add it to the seurat object (TRUE/FALSE column)
## / ! this code is not working
# scaled_matrix <- GetAssayData(seu_singlet, layer = "scale.data", assay = "SCT")
# test <- rownames(seu_singlet) %in% rownames(scaled_matrix)
# seu_singlet <- AddMetaData(seu_singlet,
#   metadata = test,
#   col.name = "sct.var.genes"
# )
# seu_singlet@meta.data$sct.var.genes

# / trying to success, with chat (no success)
# sct_vst_counts <- seu_singlet[["SCT"]]@data
# variable_features <- VariableFeatures(seu_singlet)
# variable_features
# vst_counts_variable <- sct_vst_counts[variable_features, ]
# dim(vst_counts_variable)

## / v3 problem, not relevant to SCT normalization
# seu_singlet <- readRDS("RDS/seu_singlet-clustered-27nov2024.rds")
# seu_singlet[["RNA3"]] <- as(object = seu_singlet[["RNA"]], Class = "Assay")
# DefaultAssay(seu_singlet) <- "RNA3"
# seu_singlet[["RNA"]] <- NULL
# seu_singlet <- RenameAssays(object = seu_singlet, RNA3 = "RNA")

# ! / using all genes for now, as I don't know how to exctract info for sct
seu_singlet
adata <- convertFormat(
  seu_singlet,
  from = "seurat",
  to = "anndata",
  main_layer = "counts",
  assay = "SCT",
  # assay = "ADT",
  drop_single_values = FALSE,
  outFile = "RDS/seu_singlet_sct-doublet-500-hvg_05dec24.rds.h5ad"
)
adata

# Go to python, and bring back latent.csv --------------------
latent <- read.csv("latents/latent-sct-totalVI_06dec24_wo-int-doublet-500.csv")

latent_mtx <- as.matrix(latent)
rownames(latent_mtx) <- colnames(seu_singlet)
latent_mtx <- latent_mtx[, -1] # remove 1st col
dim(latent_mtx)
colnames(latent_mtx) <- paste0("totalvi_", 1:20)
colnames(latent_mtx)

DefaultAssay(seu_singlet) <- "SCT"
seu_singlet
colnames(seu_singlet@meta.data)

seu_singlet[["totalvi.sct"]] <- CreateDimReducObject(
  embeddings = latent_mtx,
  key = "totalvi_",
  assay = DefaultAssay(seu_singlet)
)
seu_singlet
colnames(seu_singlet@meta.data)

# 9. clustering po dodaniu latent -------------------------
# / nie normalizuje bo nie trzeba i nie robie pca poniewaz używam totalvi
seu_singlet <- seu_normalize_var_scale_pca(
  seu_singlet,
  normalize = FALSE,
  cluster.name = "cluster.totalvi",
  run_pca = FALSE, # not runnign pca, changing clusters in totalvi reduction
  pca.reduction.name = "totalvi.sct", # instead of pca
  umap.reduction.name = "umap.totalvi.sct",
  dims = 1:20,
  k.param = 20,
  algorithm = 4,
  resolution = 0.1,
  group.singletons = FALSE
)

# 10. Vizualize GEX features after totalVI latent SCT clustering ---------
DefaultAssay(seu_singlet) <- "SCT"
inferno <- viridis::scale_color_viridis(option = "inferno")

vars <- t_cell_markers
vars <- macrophages
vars <- ein_like
vars <- inpc_like
vars <- c("CD8A", "CD8B")
seu_FeaturePlot(
  seu_singlet,
  assay = "SCT",
  features = vars,
  reduction = "umap.totalvi.sct",
  # reduction = "umap.sct",
  color = inferno,
  max.cutoff = NA,
  show = TRUE,
  # save_path = "plots/scv-500-doublets/inpc-like_markers",
  # ggheight = 30,
  # ggwidth = 30
)
seu_DimPlot(
  seu_singlet,
  reduction = "umap.totalvi.sct",
  # reduction = "umap.sct",
  group.by = c("MULTI_ID", "cluster.totalvi"),
  # group.by = c("MULTI_ID", "cluster.sct"),
  show = TRUE,
  # save_path = "plots/clusters",
  ggheight = NA,
  ggwidth = NA
)
seu_dittoDotPlot(
  seu_singlet,
  vars = vars,
  group.by = "seurat_clusters"
)
# 11. Vizualize ADT features -------------
DefaultAssay(seu_singlet) <- "ADT"
features <- Features(seu_singlet)
length(features)
f1 <- features[1:10]
f2 <- features[11:21]
f3 <- features[22:32]
f4 <- features[33:42]

adt_features <- c("TotalSeqC-CD8")
seu_FeaturePlot(
  seu_singlet, "ADT",
  # reduction = "umap.totalvi.sct",
  reduction = "umap.sct",
  color = inferno,
  features = adt_features,
  max.cutoff = NA,
  show = TRUE,
  # save_path = "plots/scv-500-doublets/adt-panel_4",
  ggheight = NA,
  ggwidth = NA
)

# 12. subset myeloid cluster to filter out cd8 cells -------------------
Idents(seu_singlet) <- "cluster.totalvi"
Idents(seu_singlet)

seu_singlet_myelo <- subset(seu_singlet, idents = c(2))
seu_singlet_myelo
# / 3426 samples

# / this time PCA makes more sense, bc I need to find
# / CD8-specific clusters, and totalvi latent is set
seu_singlet_myelo <- seu_normalize_var_scale_pca(
  seu_singlet_myelo,
  normalize = "sct",
  cluster.name = "cluster.pca.myelo",
  run_pca = TRUE,
  pca.reduction.name = "pca.myelo.sct",
  umap.reduction.name = "umap.myelo.pca.sct",
  dims = 1:5,
  k.param = 10,
  algorithm = 4,
  resolution = 0.4,
  group.singletons = FALSE
)
seu_DimPlot(
  seu_singlet_myelo,
  reduction = "umap.myelo.pca.sct",
  group.by = c("cluster.pca.myelo"),
  show = TRUE,
  # save_path = "plots/temp",
  ggheight = NA,
  ggwidth = NA
)
seu_FeaturePlot(
  seu_singlet_myelo,
  assay = "SCT",
  features = c("CD8A", "CD8B"),
  reduction = "umap.myelo.pca.sct",
  color = inferno,
  max.cutoff = NA,
  show = TRUE,
  # save_path = "plots/temp2",
  ggwidth = NA,
  ggheight = NA
)

saveRDS(seu_singlet, "RDS/seu_singlet_temp.rds")

# 13. remove selected clusters ---------------------
Idents(seu_singlet_myelo) <- "cluster.pca.myelo"
clusters_to_remove <- subset(seu_singlet_myelo, idents = c(6, 11))
clusters_to_remove
# / 411 samples

seu_singlet <- subset(
  seu_singlet,
  cells = setdiff(Cells(seu_singlet), Cells(clusters_to_remove))
)
seu_singlet

# 14. rereun normalize after removing CD8 cells from myelo clusters ------------
seu_singlet
seu_singlet <- seu_normalize_var_scale_pca(
  seu_singlet,
  normalize = "sct",
  cluster.name = "cluster.sct",
  run_pca = TRUE,
  pca.reduction.name = "pca.sct",
  umap.reduction.name = "umap.sct",
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

seu_DimPlot(
  seu_singlet,
  reduction = "umap.sct",
  group.by = c("cluster.sct"),
  show = TRUE,
  # save_path = "plots/temp",
  ggheight = NA,
  ggwidth = NA
)
seu_FeaturePlot(
  seu_singlet,
  assay = "SCT",
  features = c("CD8A", "CD8B"),
  reduction = "umap.sct",
  color = inferno,
  max.cutoff = NA,
  show = TRUE,
  # save_path = "plots/temp2",
  ggwidth = NA,
  ggheight = NA
)

# 15. Rerun totalVI after removing CD8 cells from myeloid cluster -----
## 15.1 convert to anndata
seu_singlet
adata <- convertFormat(
  seu_singlet,
  from = "seurat",
  to = "anndata",
  main_layer = "counts",
  # assay = "SCT",
  assay = "ADT",
  drop_single_values = FALSE,
  outFile = "RDS/seu_singlet_adt_wo-cd8inMyelo_06dec24.h5ad"
)

## 15.2 Go to python, and bring back latent.csv --------------------
latent <- read.csv("latents/latent-sct-totalVI_06dec24_after-removing-of-cd8-from-myeloid-cluster.csv")

latent_mtx <- as.matrix(latent)
rownames(latent_mtx) <- colnames(seu_singlet)
latent_mtx <- latent_mtx[, -1] # remove 1st col
dim(latent_mtx)
colnames(latent_mtx) <- paste0("totalvi_", 1:20)
colnames(latent_mtx)

DefaultAssay(seu_singlet) <- "SCT"
seu_singlet
colnames(seu_singlet@meta.data)

seu_singlet[["totalvi.sct"]] <- CreateDimReducObject(
  embeddings = latent_mtx,
  key = "totalvi_",
  assay = DefaultAssay(seu_singlet)
)
seu_singlet
colnames(seu_singlet@meta.data)

# 16. Clustering po dodaniu latent-wo-myelo -------------------------
# / nie normalizuje bo nie trzeba i nie robie pca poniewaz używam totalvi
seu_singlet <- seu_normalize_var_scale_pca(
  seu_singlet,
  normalize = FALSE,
  cluster.name = "cluster.totalvi",
  run_pca = FALSE, # not runnign pca, changing clusters in totalvi reduction
  pca.reduction.name = "totalvi.sct", # instead of pca
  umap.reduction.name = "umap.totalvi.sct",
  dims = 1:20,
  k.param = 20,
  algorithm = 4,
  resolution = 0.8,
  group.singletons = FALSE
)

# 17. Vizualize GEX features after totalVI latent SCT clustering ---------
DefaultAssay(seu_singlet) <- "SCT"
inferno <- viridis::scale_color_viridis(option = "inferno")

vars <- t_cell_markers
vars <- macrophages
vars <- ein_like
vars <- inpc_like
vars <- c("CD8A", "CD8B")
seu_FeaturePlot(
  seu_singlet,
  assay = "SCT",
  features = vars,
  reduction = "umap.totalvi.sct",
  # reduction = "umap.sct",
  color = inferno,
  max.cutoff = NA,
  show = TRUE,
  # save_path = "plots/scv-500-doublets/inpc-like_markers",
  # ggheight = 30,
  # ggwidth = 30
)
seu_DimPlot(
  seu_singlet,
  reduction = "umap.totalvi.sct",
  # reduction = "umap.sct",
  group.by = c("MULTI_ID", "cluster.totalvi"),
  # group.by = c("MULTI_ID", "cluster.sct"),
  show = TRUE,
  # save_path = "plots/clusters",
  ggheight = NA,
  ggwidth = NA
)
seu_dittoDotPlot(
  seu_singlet,
  vars = vars,
  group.by = "seurat_cluster8"
)




# 15. vizualize clonality ----------------------------------------
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
  reduction = "umap.totalvi.sct",
  group.by = "cloneSize"
) + scale_color_manual(values = rev(colorblind_vector[c(1, 3, 4, 5, 7)]))


# 16. Final Naming of clusters -------------------------------------
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

# 17. VIZ cycling cells (this should go much earlier)
test <- RunPCA(seu_singlet, features = c(s_genes, g2m_genes))
DimPlot(test, reduction = "pca")
DimPlot(test, reduction = "pca")

DimPlot(test,
  reduction = "pca",
  group.by = "Phase",
  split.by = "Phase"
)

FeaturePlot(
  test, "percent_ribo",
  reduction = "umap.sct",
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

DefaultAssay(seu_singlet) <- "SCT"
colnames(seu_singlet@meta.data)
Idents(seu_singlet) <- "cluster.totalvi"
Idents(seu_singlet)
?SetIdent
cluster_microglia.markers <- FindMarkers(
  seu_singlet,
  ident.1 = 7,
  min.pct = 0.25
)
rownames(cluster_microglia.markers)[1:100]

cluster_3.markers <- FindMarkers(
  seu_singlet,
  ident.1 = "3",
  min.pct = 0.25
)
rownames(cluster_3.markers)[1:30]

saveRDS(seu_singlet, "RDS/seu_singlet-clustered-totalvi-28nov2024.rds")
