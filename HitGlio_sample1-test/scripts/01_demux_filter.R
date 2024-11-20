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

# filter cells
seu_singlet <- PercentageFeatureSet(
  seu_singlet,
  pattern = "^MT-",
  col.name = "percent_mito"
)
seu_singlet[["log10GenesPerUmi"]] <-
  log10(seu_singlet$nFeature_RNA) / log10(seu_singlet$nCount_RNA)

Seurat::VlnPlot(
  seu_singlet,
  features = c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent_mito"
  )
)
?VlnPlot
seu_singlet <- subset(
  seu_singlet,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 7000 &
    nCount_RNA < 50000 &
    nCount_RNA > 30 &
    percent_mito < 8
)
seu_singlet
# An object of class Seurat
# 38650 features across 15250 samples within 3 assays

Seurat::VlnPlot(
  seu_singlet_subset,
  features = c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent_mito"
  )
)

# after filtering
# subset csf and tumor samples into dif. variables (for tcr analysis)
seu_singlet_csf <- subset(seu_singlet, subset = MULTI_ID == "csf")
seu_singlet_tumor <- subset(seu_singlet, subset = MULTI_ID == "tumor")
seu_singlet_csf
seu_singlet_tumor

# THIS SHOLD NOT BE DONE WITH TCR data (I think so)
# or IT MAY BE (just to vizualize clonality)
# but for gex and adt please be my guest
# analysis of tumor and csf together!!!!

# GEX analysis (both samples)
DefaultAssay(seu_singlet) <- "RNA"
seu_singlet <- NormalizeData(seu_singlet)
seu_singlet <- FindVariableFeatures(seu_singlet)
seu_singlet <- ScaleData(seu_singlet)

seu_singlet <- RunPCA(seu_singlet)
seu_singlet <- FindNeighbors(
  seu_singlet,
  reduction = "pca", dims = 1:20
)
seu_singlet <- FindClusters(
  seu_singlet,
  resolution = 0.5, verbose = FALSE
)
seu_singlet <- RunUMAP(
  seu_singlet,
  reduction = "pca", dims = 1:20
)

DimPlot(
  seu_singlet,
  group.by = "MULTI_ID", reduction = "umap"
)

colnames(seu_singlet@meta.data)
DimPlot(
  seu_singlet,
  group.by = "seurat_clusters", reduction = "umap",
  label = TRUE
)

library(viridis)
mycolor2 <- scale_color_viridis(option = "viridis")

FeaturePlot(
  seu_singlet, c("CD8B", "NKG7"),
  # max.cutoff = 2,
  label = TRUE
) & mycolor2

FeaturePlot(
  seu_singlet, c("CD4"),
  # max.cutoff = 2,
  label = TRUE
) & mycolor2

FeaturePlot(
  seu_singlet, c("CX3CR1"),
  # max.cutoff = 2
) & mycolor2

FeaturePlot(
  seu_singlet, c("H3F3A", "IDH1", "IDH2", "SETD2"),
  # max.cutoff = 2,
  label = TRUE
) & mycolor2

# cell metaprograms
# from: Liu, Ilon, Li Jiang, Erik R. Samuelsson, Sergio Marco Salas, Alexander Beck, Olivia A. Hack, Daeun Jeong, et al. “The Landscape of Tumor Cell States and Spatial Organization in H3-K27M Mutant Diffuse Midline Glioma across Age and Location.” Nature Genetics 54, no. 12 (December 2022): 1881–94. https://doi.org/10.1038/s41588-022-01236-3. # nolint
ac_like <- c("CLU", "AQP4", "AGT", "SPARCL1", "VIM", "GFAP", "CRYAB", "MLC1", "PLTP", "LAMB2", "CD99", "APOE", "S1PR1", "SPARC", "GJA1", "EDNRB", "SLC1A3", "F3", "CD44", "TNC", "C4A", "EZR", "HEPN1", "ID3", "CCDC80", "S100A10", "IFI6", "RASSF4", "ATP1A2", "ALDOC") # nolint
mes_like <- c("VIM", "S100A10", "TIMP1", "GAP43", "CLU", "S100A6", "TAGLN2", "TMSB10", "MYL12A", "S100A16", "DDIT3", "SPP1", "LMNA", "TNFRSF12A", "VMP1", "IL13RA2", "SQSTM1", "SLC2A3", "OCIAD2", "CD63", "SPARCL1", "LGALS1", "MT2A", "SCG2", "S100A11", "CAV1", "CD59", "EMP1") # nolint
# mes like not available in out dataset
# "RP11-430H10.1"
# "ORAOV1"

# cell types from: Neftel, Cyril, Julie Laffy, Mariella G. Filbin, Toshiro Hara, Marni E. Shore, Gilbert J. Rahme, Alyssa R. Richman, et al. “An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma.” Cell 178, no. 4 (August 8, 2019): 835-849.e21. https://doi.org/10.1016/j.cell.2019.06.024. #nolint
oligoden <- c("MBP", "TF", "PLP1", "MAG", "MOG", "CLDN11")
tcels <- c("CD2", "CD3D", "CD3E", "CD3G")
macrophages <- c("CD14", "AIF1", "FCER1G", "FCGR3A", "TYROBP", "CSF1R")

fp1 <- FeaturePlot(
  seu_singlet, mes_like,
  # max.cutoff = 2,
  label = TRUE
) & mycolor2

library(dittoSeq)
dittoSeq::dittoDotPlot(
  seu_singlet,
  vars = mes_like,
  group.by = "MULTI_ID"
  # group.by = "seurat_clusters"
)



library(plot1cell)


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
