# according to: https://satijalab.org/seurat/articles/hashing_vignette

library(Seurat)
library(ggplot2)
library(Matrix)
library(tidyverse)
library(clustree)

# setwd("/mnt/sda4/singleCell_LAB/HitGlio_sample1-test")

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
joint_bcs

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

# seu <- HTODemux(seu, assay = "HTO", positive.quantile = 0.99)
# table(seu$HTO_classification.global)

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

seu_demux

seu_singlet <- subset(seu_demux, idents = c("tumor", "csf"))
head(seu_singlet@meta.data)

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

DefaultAssay(seu_joined_singlet) <- "RNA"
p2 <- FeaturePlot(seu_joined_singlet, "Itga4", max.cutof = 2, cols = mycolor)

p1 | p2


library(ProjecTILs)
library(scRepertoire)
library(tidyverse)

tcr_sample <- read_csv("/mnt/sdb1/runs/sample1_multilane_spec_index_NNN/outs/per_sample_outs/sample1/vdj_t/filtered_contig_annotations.csv")

head(tcr_sample)
dim(tcr_sample)

tcr_sample_list <- list(tcr_sample)
tcr_sample_list

combined_tcr <- combineTCR(
  tcr_sample_list,
  # to make things easier, w/o 'samples' argument
  # it wont add samples names to barcodes
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

seu_singlet_tcr <- combineExpression(
  combined_tcr,
  seu_singlet,
  cloneCall = "gene",
  proportion = TRUE
)
colnames(seu_singlet_tcr@meta.data)

colorblind_vector <- hcl.colors(n = 7, palette = "inferno", fixup = TRUE)

DimPlot(
  seu_singlet_tcr,
  group.by = "cloneSize"
) + scale_color_manual(values = rev(colorblind_vector[c(1, 3, 4, 5, 7)]))


# CD4
ref_cd4 <- load.reference.map(
  "singleCell_LAB/scREP_projectTIL_tut/data/refs/pTILs_hsa_cd4t.rds"
)

seu_singlet_tcr <- Run.ProjecTILs(
  seu_singlet_tcr,
  ref = ref_cd4,
  ncores = 1
)

p1 <- plot.projection(ref_cd4) + theme(aspect.ratio = 1)
p2 <- plot.projection(
  ref_cd4, seu_singlet_tcr,
  linesize = 0.3, pointsize = 0.5
) + theme(aspect.ratio = 1)
p1 | p2

plot.statepred.composition(ref_cd4, seu_singlet_tcr, metric = "Percent")

plot.states.radar(
  ref = ref_cd4, seu_singlet_tcr, min.cells = 30
)

colnames(seu_singlet_tcr@meta.data)

# CD8
ref_cd8 <- load.reference.map(
  "singleCell_LAB/scREP_projectTIL_tut/data/refs/pTILs_hsa_cd8t.rds"
)

seu_tcr <- Run.ProjecTILs(
  seu_singlet_tcr,
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

seu_tcr
