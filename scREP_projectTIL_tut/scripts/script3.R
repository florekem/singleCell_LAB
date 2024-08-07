# https://carmonalab.github.io/ProjecTILs_CaseStudies/Xiong19_TCR.html

options(timeout = 5000)
options(warn = 1)

library(renv)
renv::restore()

library(Seurat)
library(ggplot2)
library(patchwork)
library(ProjecTILs)
library(scRepertoire)

setwd("/mnt/sda4/singleCell_LAB/scREP_projectTIL_tut/data/MC38_TILs/")

# --- prepare gex data ----------------------------------------------------
sparse_matrix <- Read10X("data/gex")

seu <- CreateSeuratObject(
  sparse_matrix,
  project = "Xiong_TIL",
  min.cells = 3,
  min.features = 50
)
seu
# An object of class Seurat
# 21075 features across 6182 samples within 1 assay
# Active assay: RNA (21075 features, 0 variable features)
#  1 layer present: counts
head(seu@meta.data)
rownames(seu@meta.data) # == colnames(seu)
# [6178] "TTTGGTTCAATGAAAC-7" "TTTGGTTCATACAGCT-7" "TTTGGTTTCTTAGAGC-7"
# [6181] "TTTGTCAAGTTGAGAT-7" "TTTGTCACAGGATTGG-7"
seu$sampleID <- substring(rownames(seu@meta.data), 18)
table(seu@meta.data$sampleID)
#   4    5    6    7
# 1759 1746 1291 1386
sample_ids <- c(4:7)
sample_names <- c("Mouse_1", "Mouse_2", "Mouse_3", "Mouse_4")
names(sample_names) <- sample_ids
sample_names
#         4         5         6         7
# "Mouse_1" "Mouse_2" "Mouse_3" "Mouse_4"
seu$sampleName <- factor(
  seu$sampleID,
  levels = sample_ids, labels = sample_names # smart af
)
table(seu@meta.data$sampleName)
# Mouse_1 Mouse_2 Mouse_3 Mouse_4
#    1759    1746    1291    1386

# --- prepare tcr data ----------------------------------------------------
