# https://carmonalab.github.io/ProjecTILs_CaseStudies/Xiong19_TCR.html

options(timeout = 5000)
options(warn = 1)

library(Seurat)
library(ggplot2)
library(patchwork)
library(ProjecTILs)
library(scRepertoire)
library(tidyverse)

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
seu$sample <- substring(rownames(seu@meta.data), 18)
table(seu@meta.data$sample)
#   4    5    6    7
# 1759 1746 1291 1386
sample_ids <- c(4:7)
sample_names <- c("Mouse_1", "Mouse_2", "Mouse_3", "Mouse_4")
names(sample_names) <- sample_ids
sample_names
#         4         5         6         7
# "Mouse_1" "Mouse_2" "Mouse_3" "Mouse_4"
seu$sampleName <- factor(
  seu$sample,
  levels = sample_ids, labels = sample_names # smart af
)
table(seu@meta.data$sampleName)
# Mouse_1 Mouse_2 Mouse_3 Mouse_4
#    1759    1746    1291    1386

# resolve barcodes (to match never version of scRep)
# tutorial is outdated
seu <- RenameCells(seu, seu@meta.data$sample)
head(rownames(seu@meta.data))
# --- prepare tcr data ----------------------------------------------------
# make list with tcr samples
sample_ids_vdj <- c(4:7)
names(sample_ids_vdj) <- c(35:38)
sample_ids_vdj
# 35 36 37 38
#  4  5  6  7
names(sample_ids_vdj)
# [1] "35" "36" "37" "38"
vdj_samples_list <- list()
# miejsce na liście musi zaczynać się od 1 i kończyć na 4, an nie od 35,
# dlatego trzeba taki myk z names(sample_ids_vdj)[i]
for (i in seq_along(sample_ids_vdj)) { # i = 1,2,3,4
  print(i)
  sample_name <- names(sample_ids_vdj)[i] # 1 == 35, 2 = 36 etc...
  vdj_samples_list[[i]] <- read.csv(
    sprintf(
      "data/tcr/filtered_contig_annotations_%s.csv", sample_name
    )
  )
  # resolving the problem with barcodes (as mentioned below)
  vdj_samples_list[[i]]$barcode <- sub(
    "\\d$", sample_ids_vdj[[i]], vdj_samples_list[[i]]$barcode
    # will I ever learn regex?
    # \d == [0-9]; $ == end of line
    # additional "\" is an escape symbol in R
  )
  # adding sampleID to "raw_clonotype_id" (mentioned below)
  # how important is that?
  vdj_samples_list[[i]]$raw_clonotype_id <- paste0(
    vdj_samples_list[[i]]$raw_clonotype_id,
    "-",
    sample_ids_vdj[[i]]
  )
}
# the problem with barcodes:
head(vdj_samples_list[[1]])
# AAACCTGAGGCAATTA-1
head(vdj_samples_list[[2]])
# AAACCTGAGCAATCTC-1
# meanwhile...
head(seu@meta.data)
# AAACGGGCATCCTAGA-4
tail(seu@meta.data)
# TTTGCGCCAGCATACT-7
# the problem with "raw_clonotype_id"
# clonotype194 etc... zmianiamy na clonotype194-5 etc...

# --- combine chains TCR --------------------------------------------------
combined <- combineTCR(
  vdj_samples_list,
  samples = sample_ids_vdj,
  # ID = names(sample_ids_vdj),
  # cells = "T-AB", # ??
  removeNA = TRUE,
  removeMulti = TRUE
)
head(combined[[4]])
colnames(combined[[1]])
# barcode:
# 4_35_AAACCTGTCGCGCCAA-4
# barcode bez "samples" i "ID" w combineTCR
# AAACCTGTCGCGCCAA-4
# może można bez tego jednak (jeśli jedyne co to robi to dodaje cyferki?)
seu <- combineExpression(
  combined,
  seu,
  cloneCall = "gene",
  group.by = "sample",
  proportion = TRUE,
)
colnames(seu@meta.data)
head(seu@meta.data)
# cloneSize = c(
#   Single = 1, Small = 5, Medium = 20, Large = 100,
#   Hyperexpanded = 500
# )
refer <- load.reference.map()
refer
seu <- NormalizeData(seu)
seu <- Run.ProjecTILs(
  seu,
  ref = refer,
  ncores = 4
)
seu@meta.data

BPPARAM <- BiocParallel::MulticoreParam(workers = 4)
register(SerialParam())
library(BiocParallel)
