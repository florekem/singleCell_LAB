library(ProjecTILs)
library(scRepertoire)
library(tidyverse)

setwd("/mnt/sda4/singleCell_LAB/HitGlio_sample1-test")

tcr_sample <- read_csv("/mnt/sdb1/runs/sample1_multilane_spec_index_NNN/outs/per_sample_outs/sample1/vdj_t/filtered_contig_annotations.csv")


# analysis of tcr only for tumor subset
DefaultAssay(seu_singlet_tumor) <- "RNA"
seu_singlet_tumor <- NormalizeData(seu_singlet_tumor)
seu_singlet_tumor <- FindVariableFeatures(seu_singlet_tumor)
seu_singlet_tumor <- ScaleData(seu_singlet_tumor)

seu_singlet_tumor <- RunPCA(seu_singlet_tumor)
seu_singlet_tumor <- FindNeighbors(
  seu_singlet_tumor,
  reduction = "pca", dims = 1:50
)
seu_singlet_tumor <- FindClusters(
  seu_singlet_tumor,
  resolution = 1.1, verbose = FALSE
)
seu_singlet_tumor <- RunUMAP(
  seu_singlet_tumor,
  reduction = "pca", dims = 1:30
)

# combine tcr for csf only
# keep only barcodes precest in csf sample
tcr_sample_tumor <- tcr_sample[
  tcr_sample$barcode %in% colnames(seu_singlet_tumor),
]
dim(tcr_sample_tumor)
# 13022 31

tcr_sample_tumor_list <- list(tcr_sample_tumor)

combined_tumor_tcr <- combineTCR(
  tcr_sample_tumor_list,
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE
)

seu_singlet_tumor <- combineExpression(
  combined_tumor_tcr,
  seu_singlet_tumor,
  cloneCall = "gene",
  proportion = TRUE
)

# CD4 tumor only
ref_cd4 <- load.reference.map(
  "/mnt/sda4/singleCell_LAB/scREP_projectTIL_tut/data/refs/pTILs_hsa_cd4t.rds"
)

seu_tumor_tcr4 <- Run.ProjecTILs(
  seu_singlet_tumor,
  ref = ref_cd4,
  ncores = 1
)
seu_tumor_tcr4

p1 <- plot.projection(ref_cd4) + theme(aspect.ratio = 1)
p2 <- plot.projection(
  ref_cd4, seu_tumor_tcr4,
  linesize = 0.3, pointsize = 0.5
) + theme(aspect.ratio = 1)
p1 | p2

ggsave("plots/plot1.png", p1 | p2)

plot.statepred.composition(ref_cd4, seu_tumor_tcr4, metric = "Percent")

plot.states.radar(
  ref = ref_cd4, seu_tumor_tcr4, min.cells = 30
)

# CD8 tumor only
ref_cd8 <- load.reference.map(
  "/mnt/sda4/singleCell_LAB/scREP_projectTIL_tut/data/refs/pTILs_hsa_cd8t.rds"
)

seu_tumor_tcr8 <- Run.ProjecTILs(
  seu_singlet_tumor,
  ref = ref_cd8,
  ncores = 1
)
seu_tumor_tcr8

p1 <- plot.projection(ref_cd8) + theme(aspect.ratio = 1)
p2 <- plot.projection(
  ref_cd8, seu_tumor_tcr8,
  linesize = 0.3, pointsize = 0.5
) + theme(aspect.ratio = 1)
p1 | p2

ggsave("plots/plot1.png", p1 | p2)

plot.statepred.composition(ref_cd8, seu_tumor_tcr8, metric = "Percent")

plot.states.radar(
  ref = ref_cd8, seu_tumor_tcr8, min.cells = 30
)
