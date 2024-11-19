library(ProjecTILs)
library(scRepertoire)
library(tidyverse)


setwd("/mnt/sda4/singleCell_LAB/HitGlio_sample1-test")

tcr_sample <- read_csv("/mnt/sdb1/runs/sample1_multilane_spec_index_NNN/outs/per_sample_outs/sample1/vdj_t/filtered_contig_annotations.csv")

head(tcr_sample)
dim(tcr_sample)

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


# CD4 together
ref_cd4 <- load.reference.map(
  "singleCell_LAB/scREP_projectTIL_tut/data/refs/pTILs_hsa_cd4t.rds"
)

seu_tcr4 <- Run.ProjecTILs(
  seu_singlet_tcr,
  ref = ref_cd4,
  ncores = 1
)

p1 <- plot.projection(ref_cd4) + theme(aspect.ratio = 1)
p2 <- plot.projection(
  ref_cd4, seu_tcr4,
  linesize = 0.3, pointsize = 0.5
) + theme(aspect.ratio = 1)
p1 | p2

ggsave("plots/plot3.png", p1 | p2)

plot.statepred.composition(ref_cd4, seu_tcr4, metric = "Percent")

plot.states.radar(
  ref = ref_cd4, seu_tcr4, min.cells = 30
)

colnames(seu_singlet_tcr@meta.data)

# CD8 together
ref_cd8 <- load.reference.map(
  "singleCell_LAB/scREP_projectTIL_tut/data/refs/pTILs_hsa_cd8t.rds"
)

seu_tcr8 <- Run.ProjecTILs(
  seu_singlet_tcr,
  ref = ref_cd8,
  ncores = 1
)

p1 <- plot.projection(ref_cd8) + theme(aspect.ratio = 1)
p2 <- plot.projection(
  ref_cd8, seu_tcr8,
  linesize = 0.3, pointsize = 0.5
) + theme(aspect.ratio = 1)
p1 | p2

plot.statepred.composition(ref_cd8, seu_tcr8, metric = "Percent")

plot.states.radar(
  ref = ref_cd8, seu_tcr8, min.cells = 30
)
