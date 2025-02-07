library(ProjecTILs)
library(scRepertoire)
library(tidyverse)

setwd("/mnt/sda4/singleCell_LAB/HitGlio_sample1-test")

tcr_sample <- read_csv("/mnt/sdb1/runs/sample1_multilane_spec_index_NNN/outs/per_sample_outs/sample1/vdj_t/filtered_contig_annotations.csv")


# I should create variable here to indicate if its csf or csf
# for now just use:
# :%s/csf/csf/g
# :%s/csf/csf/g

#
# analysis of tcr only for csf subset
DefaultAssay(seu_singlet_csf) <- "RNA"
seu_singlet_csf <- NormalizeData(seu_singlet_csf)
seu_singlet_csf <- FindVariableFeatures(seu_singlet_csf)
seu_singlet_csf <- ScaleData(seu_singlet_csf)
seu_singlet_csf <- RunPCA(seu_singlet_csf)
seu_singlet_csf <- FindNeighbors(
  seu_singlet_csf,
  reduction = "pca", dims = 1:20
)
seu_singlet_csf <- FindClusters(
  seu_singlet_csf,
  resolution = 1.1, verbose = FALSE
)
seu_singlet_csf <- RunUMAP(
  seu_singlet_csf,
  reduction = "pca", dims = 1:20
)

# combine tcr for csf only
# keep only barcodes precest in csf sample
tcr_sample_csf <- tcr_sample[
  tcr_sample$barcode %in% colnames(seu_singlet_csf),
]
dim(tcr_sample_csf)
# 2218 csf
# 12687 csf
# 2225 31 było bez filtrowania. czy to nie za duzo?

tcr_sample_csf_list <- list(tcr_sample_csf)

combined_csf_tcr <- combineTCR(
  tcr_sample_csf_list,
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE
)

seu_singlet_csf <- combineExpression(
  combined_csf_tcr,
  seu_singlet_csf,
  cloneCall = "gene",
  proportion = TRUE
)

# CD4
ref_cd4 <- load.reference.map(
  "/mnt/sda4/singleCell_LAB/scREP_projectTIL_tut/data/refs/pTILs_hsa_cd4t.rds"
)

seu_csf_tcr4 <- Run.ProjecTILs(
  seu_singlet_csf,
  ref = ref_cd4,
  ncores = 1
)
seu_csf_tcr4

p1 <- plot.projection(ref_cd4) + theme(aspect.ratio = 1)
p2 <- plot.projection(
  ref_cd4, seu_csf_tcr4,
  linesize = 0.3, pointsize = 0.5
) + theme(aspect.ratio = 1)
p1 | p2

ggsave("plots/seu_csf_tcr4.png", p1 | p2)

p3 <- plot.statepred.composition(ref_cd4, seu_csf_tcr4, metric = "Percent")
p3 + theme(
  legend.text = element_text(size = 15),
  plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 14)
) + labs(
  title = "CD4 scf",
  x = "Cell state",
  y = "Percentage of cells"
) + scale_y_continuous(limits = c(0, 60))




ggsave("plots/seu_csf_tcr4_percent-composition.png", p3)

p4 <- plot.states.radar(
  ref = ref_cd4, seu_csf_tcr4, min.cells = 30
)
ggsave("plots/seu_csf_tcr4_genes-radar.png", p4)

ref_cd8 <- load.reference.map(
  "/mnt/sda4/singleCell_LAB/scREP_projectTIL_tut/data/refs/pTILs_hsa_cd8t.rds"
)


seu_csf_tcr8 <- Run.ProjecTILs(
  seu_singlet_csf,
  ref = ref_cd8,
  ncores = 1
)
seu_csf_tcr8

p5 <- plot.projection(ref_cd8) + theme(aspect.ratio = 1)
p6 <- plot.projection(
  ref_cd8, seu_csf_tcr8,
  linesize = 0.3, pointsize = 0.5
) + theme(aspect.ratio = 1)
p5 | p6

ggsave("plots/seu_csf_tcr8.png", p5 | p6)

p7 <- plot.statepred.composition(ref_cd8, seu_csf_tcr8, metric = "Percent")
p7 + theme(
  legend.text = element_text(size = 15),
  plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 14)
) + labs(
  title = "CD8 CSF",
  x = "Cell state",
  y = "Percentage of cells"
)


ggsave("plots/seu_csf_tcr8_percent-composition.png", p7)

p8 <- plot.states.radar(
  ref = ref_cd8, seu_csf_tcr8, min.cells = 30
)
ggsave("plots/seu_csf_tcr8_genes-radar.png", p8)
