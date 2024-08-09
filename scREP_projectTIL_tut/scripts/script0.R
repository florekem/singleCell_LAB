library(patchwork)
library(ggplot2)
library(reshape2)
library(Seurat)
library(ProjecTILsu)


options(timeout = 5000)
ddir <- "input/Bassez_breast_data"
ddir <- "/mnt/sda4/singleCell_LAB/scREP_projectTIL_tut/scripts/input/Bassez_breast_data"

dir.create(ddir)
dataUrl <- "https://figshare.com/ndownloader/files/47900725"
download.file(dataUrl, paste0(ddir, "/tmp.zip"))
unzip(paste0(ddir, "/tmp.zip"), exdir = ddir)
file.remove(paste0(ddir, "/tmp.zip"))

# Count matrices
f1 <- sprintf("%s/1864-counts_tcell_cohort1.rds", ddir)

cohort1 <- readRDS(f1)
dim(cohort1)

meta1 <- read.csv(sprintf("%s/1870-BIOKEY_metaData_tcells_cohort1_web.csv", ddir))
rownames(meta1) <- meta1$Cell

data.seurat <- CreateSeuratObject(cohort1, project = "Cohort1_IT", meta.data = meta1)
data.seurat <- NormalizeData(data.seurat)

data.seurat <- subset(data.seurat, subset = timepoint == "Pre")

data.seurat <- subset(data.seurat, subset = cellSubType %in% c(
    "NK_REST", "Vg9Vd2_gdT",
    "gdT", "NK_CYTO"
), invert = T)

ds <- 1000
min.cells <- 200

tab <- table(data.seurat$patient_id)
keep <- names(tab)[tab > min.cells]
data.seurat <- subset(data.seurat, patient_id %in% keep)


Idents(data.seurat) <- "patient_id"
data.seurat <- subset(data.seurat, cells = WhichCells(data.seurat, downsample = ds))
table(data.seurat$patient_id)

download.file("https://figshare.com/ndownloader/files/41414556", destfile = "CD8T_human_ref_v1.rds")

ref.cd8 <- load.reference.map("CD8T_human_ref_v1.rds")
ref.cd8 <- readRDS("CD8T_human_ref_v1.rds")
ref.cd8 <- UpdateSeuratObject(ref.cd8)
ncores <- 8

BPPARAM = BiocParallel::MulticoreParam(workers = 1)

DefaultAssay(ref.cd8) <- "integrated"
data.seurat <- ProjecTILs.classifier(data.seurat, ref.cd8, ncores = 1, split.by = "patient_id")

head(data.seurat@meta.data)
