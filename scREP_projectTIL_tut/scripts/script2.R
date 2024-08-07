# https://carmonalab.github.io/ProjecTILs_CaseStudies/Bassez_BC.html

library(patchwork)
library(ggplot2)
library(reshape2)
library(Seurat)
library(ProjecTILs)

options(timeout = 5000)

setwd("/mnt/sda4/singleCell_LAB/scREP_projectTIL_tut/")
# --- get data -----------------------------------------------------------------
ddir <- "data/Bassez_breast_data"
if (!file.exists(ddir)) {
  dir.create(ddir)
  data_url <- "https://figshare.com/ndownloader/files/47900725"
  download.file(data_url, paste0(ddir, "/tmp.zip"))
  unzip(paste0(ddir, "/tmp.zip"), exdir = ddir)
  file.remove(paste0(ddir, "/tmp.zip"))
}

# --- Count matrices -----------------------------------------------------------
f1 <- sprintf("%s/1864-counts_tcell_cohort1.rds", ddir)
cohort1 <- readRDS(f1) # sparse matrix
dim(cohort1)

# --- Metadata -----------------------------------------------------------------
meta1 <- read.csv(
  sprintf("%s/1870-BIOKEY_metaData_tcells_cohort1_web.csv", ddir)
)
rownames(meta1) <- meta1$Cell
head(meta1)
# --- seurat object -----------------------------------------------------------
seurat <- CreateSeuratObject(
  cohort1,
  project = "Cohort1_IT", meta.data = meta1
)
seurat <- NormalizeData(seurat)
seurat
# An object of class Seurat
# 25288 features across 53382 samples within 1 assay
# Active assay: RNA (25288 features, 0 variable features)
#  2 layers present: counts, data
colnames(seurat@meta.data)
#  [1] "orig.ident"   "nCount_RNA"   "nFeature_RNA" "Cell"  "patient_id"
#  [6] "timepoint"    "expansion"    "BC_type"  "cellSubType"  "cohort"

# --- Subset pre-treatment biopsies ------------------------------------------
table(seurat$timepoint)
#  On   Pre
# 27528 25854
seurat_pre <- subset(
  seurat,
  subset = timepoint == "Pre"
)

# --- Subset T cells according to annotation by the authors -----------------
table(seurat_pre$cellSubType)
#                CD4_EM                CD4_EX  CD4_EX_Proliferating
#                  4810                  1707                   107
#                 CD4_N               CD4_REG CD4_REG_Proliferating
#                  3007                  2861                    69
#                CD8_EM              CD8_EMRA                CD8_EX
#                  5443                   290                  2286
#  CD8_EX_Proliferating                 CD8_N                CD8_RM
#                   340                   354                  1963
#                   gdT               NK_CYTO               NK_REST
#                  1000                   230                   938
#            Vg9Vd2_gdT
#                   449
seurat_pre_t <- subset(
  seurat_pre,
  subset = cellSubType %in% c(
    "NK_REST", "Vg9Vd2_gdT",
    "gdT", "NK_CYTO"
  ), invert = TRUE # clever
)
head(seurat_pre_t@meta.data)

# --- Downsampling ---------------------------------------------------------
# We will remove samples that are too small and downsample large ones,
# in order to have similar contribution from all patients/samples.
# Downsampling large samples also speeds up downstream computations.
table(seurat_pre_t$patient_id)
#  BIOKEY_1 BIOKEY_10 BIOKEY_11 BIOKEY_12 BIOKEY_13 BIOKEY_14 BIOKEY_15 BIOKEY_16 #nolint
#      3062      1008       555      2531      1294       939      1036      2120 #nolint
# BIOKEY_17 BIOKEY_18 BIOKEY_19  BIOKEY_2 BIOKEY_20 BIOKEY_21 BIOKEY_22 BIOKEY_23 #nolint
#        97       196      1440       960       110       203        68        35 #nolint
# BIOKEY_24 BIOKEY_25 BIOKEY_26 BIOKEY_27 BIOKEY_28 BIOKEY_29  BIOKEY_3 BIOKEY_30 #nolint
#       347        99        86       601       594       163       266        61 #nolint
# BIOKEY_31  BIOKEY_4  BIOKEY_5  BIOKEY_6  BIOKEY_7  BIOKEY_8  #nolint
#       349      2326      1700       616        63       140
min.cells <- 200
tab <- table(seurat_pre_t$patient_id) # remove patients with low number of cells
keep <- names(tab)[tab > min.cells]
seurat_pre_t <- subset(
  seurat_pre_t, patient_id %in% keep # need to remember that (keep == chr[1:19])
)

ds <- 1000
Idents(seurat_pre_t) <- "patient_id"
seurat_pre_t <- subset(
  seurat_pre_t,
  cells = WhichCells(seurat_pre_t, downsample = ds) # clever
)
table(seurat_pre_t$patient_id)
#  BIOKEY_1 BIOKEY_10 BIOKEY_11 BIOKEY_12 BIOKEY_13 BIOKEY_14 BIOKEY_15 BIOKEY_16  #nolint
#      1000      1000       555      1000      1000       939      1000      1000  #nolint
# BIOKEY_19  BIOKEY_2 BIOKEY_21 BIOKEY_24 BIOKEY_27 BIOKEY_28  BIOKEY_3 BIOKEY_31  #nolint
#      1000       960       203       347       601       594       266       349  #nolint
#  BIOKEY_4  BIOKEY_5  BIOKEY_6  #nolint
#      1000      1000       616  #nolint

# --- ProjecTILs analysis ----------------------------------------------------
# --- Load reference maps as seurat objects ----------------------------------
# CD8
download.file(
  "https://figshare.com/ndownloader/files/41414556",
  destfile = "CD8T_human_ref_v1.rds" # seurat object
)
ref_cd8 <- load.reference.map("scripts/CD8T_human_ref_v1.rds")
class(ref_cd8)
# [1] "SeuratObject"
ref_cd8
# An object of class Seurat
# 55614 features across 10045 samples within 3 assays
# Active assay: RNA (27407 features, 0 variable features)
#  2 layers present: counts, data
#  2 other assays present: integrated, RNA_trimmed
#  2 dimensional reductions calculated: pca, umap
colnames(ref_cd8@meta.data)
dim(ref_cd8)
head(ref_cd8)

# CD4
download.file(
  "https://figshare.com/ndownloader/files/39012395",
  destfile = "CD4T_human_ref_v1.rds"
)
ref_cd4 <- load.reference.map("CD4T_human_ref_v1.rds")

a <- DimPlot(
  ref_cd8,
  cols = ref_cd8@misc$atlas.palette, label = TRUE
) + theme(aspect.ratio = 1) +
  ggtitle("CD8 T reference") + NoLegend()

b <- DimPlot(
  ref_cd4,
  cols = ref_cd4@misc$atlas.palette, label = TRUE
) + theme(aspect.ratio = 1) +
  ggtitle("CD4 T reference") + NoLegend()

a | b

# --- Show key marker genes for the reference subtypes ------------------------
DefaultAssay(ref_cd8) <- "RNA"
Idents(ref_cd8) <- "functional.cluster"

genes <- c(
  "SELL", "TCF7", "LEF1", "CCR7", "S1PR1", "LMNA", "IL7R", "GZMK", "FGFBP2",
  "FCGR3A", "XCL1", "XCL2", "CD200", "CRTAM", "GNG4", "TOX", "PDCD1", "HAVCR2",
  "GZMB", "PRF1", "LAG3", "KLRB1", "TRAV1-2"
)

DotPlot(
  ref_cd8,
  features = genes, cols = "RdBu", scale = TRUE, col.max = 1.5
) + theme(axis.text.x = element_text( # need to remember that!
  angle = 45,
  hjust = 1
)) + ggtitle("CD8 T cell reference map")

DefaultAssay(ref_cd4) <- "RNA"
Idents(ref_cd4) <- "functional.cluster"

genes <- c(
  "SELL", "TCF7", "S1PR1", "KLF2", "IL7R", "CCL5", "GZMK", "EOMES", "TBX21",
  "CXCR6", "GZMA", "GZMB", "PRF1", "GNLY", "CXCL13", "PDCD1", "TOX", "LAG3",
  "HAVCR2", "KLRB1", "IL17A", "FOXP3", "IL2RA"
)

DotPlot(
  ref_cd4,
  features = genes, cols = "RdBu", scale = TRUE, col.max = 1.5
) + theme(axis.text.x = element_text(
  angle = 45,
  hjust = 1
)) + ggtitle("CD4 T cell reference map")

# --- Reference-based analysis -----------------------------------------------
# Classify CD8 T subtypes
ncores <- 12

DefaultAssay(ref_cd8) <- "integrated"
DefaultAssay(seurat_pre_t)
# [1] "RNA"
# [I would definitely recommend NOT to use the integrated assay for your query data.] #nolint
# https://github.com/carmonalab/ProjecTILs/issues/72

seurat_pre_t <- ProjecTILs.classifier(
  seurat_pre_t, ref_cd8,
  split.by = "patient_id" # For optimal batch-effect correction, we recommend projecting each patient/batch separately (split.by) #nolint
)
table(seurat_pre_t$functional.cluster, useNA = "ifany")

ref <- load.reference.map()
ref
ref_cd8
