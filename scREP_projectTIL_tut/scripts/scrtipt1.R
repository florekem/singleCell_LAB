library(scRepertoire)
library(ProjecTILs)
library(Seurat)
# ----------------------------------------------------------------------
# The built-in example data is derived from the 10x Cell Ranger
# pipeline, so it is ready to go for
# downstream processing and analysis.
data("contig_list") # the data built into scRepertoire
head(contig_list[[1]]) # list of lists
colnames(contig_list[[1]])
#  [1] "barcode"   "is_cell"   "contig_id"   "high_confidence"
#  [5] "length"      "chain"            "v_gene"           "d_gene"
#  [9] "j_gene"   "c_gene"    "full_length"  "prooductive"
# [13] "cdr3" cdr3_nt" reads"   "umis"
# [17] "raw_clonotype_id" "raw_consensus_id"
dim(contig_list[[1]])
# [1] 5504   18
# ----------------------------------------------------------------------
# combining contings into clones:
# The first step in getting clones is
# to use the single-cell barcodes to organize cells into paired sequences.
combined_tcr <- combineTCR(
  contig_list,
  samples = c(
    "P17B", "P17L", "P18B", "P18L",
    "P19B", "P19L", "P20B", "P20L"
  ),
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE
)
class(combined_tcr) # list
class(combined_tcr[1]) # list
class(combined_tcr[[1]]) # data frame
str(combined_tcr)
head(combined_tcr[[1]])
colnames(combined_tcr[[1]])
# [1] "barcode"  "sample"   "TCR1" "cdr3_aa1" "cdr3_nt1" "TCR2"
# [7] "cdr3_aa2" "cdr3_nt2" "CTgene" "CTnt" "CTaa"  "CTstrict"
dim(combined_tcr[[1]])
# [1] 2805   13
# ----------------------------------------------------------------------
# (optional) addVariable
combined_tcr <- addVariable(
  combined_tcr,
  variable.name = "Type",
  variables = rep(c("B", "L"), 4)
)
colnames(combined_tcr[[1]])
# [1] "barcode"  "sample"   "TCR1"     "cdr3_aa1" "cdr3_nt1" "TCR2"
# [7] "cdr3_aa2" "cdr3_nt2" "CTgene"   "CTnt"     "CTaa"     "CTstrict"
# [13] "Type"
head(rownames(combined_tcr[[1]]))
# [1] "1"  "3"  "5"  "7"  "9"  "10"
# ----------------------------------------------------------------------
# (optional) subsetClones
subset1 <- subsetClones(
  combined_tcr,
  name = "sample",
  variables = c("P18L", "P18B")
)
head(subset1[[1]])
# Alternatively, we can also just select the list elements
# after combineTCR() or combineBCR().
subset2 <- combined_tcr[c(3, 4)]
head(subset2[[1]])
# exportClones (to save for later use or to use in other pipelines
exportClones(
  combined,
  write.file = TRUE,
  dir = "~/Documents/MyExperiment/Sample1/",
  file.name = "clones.csv"
)
# Basic Clonal Visualizations
clonalQuant(
  combined_tcr,
  cloneCall = "strict",
  chain = "both",
  scale = TRUE
)
clonalQuant(
  combined_tcr,
  cloneCall = "gene", group.by = "Type", scale = TRUE
)
# ----------------------------------------------------------------------
# Combining Clones and Single-Cell Objects
data(scRep_example)
scRep_example
# An object of class Seurat
# 2000 features across 500 samples within 1 assay
# Active assay: RNA (2000 features, 2000 variable features)
#  2 layers present: counts, data
#  1 dimensional reduction calculated: umap
seu_tcr <- combineExpression(
  combined_tcr,
  scRep_example, # seu object
  cloneCall = "gene",
  group.by = "sample",
  proportion = TRUE
)
colnames(seu_tcr@meta.data)
#  [1] "orig.ident"       "nCount_RNA"       "nFeature_RNA"     "mito.genes"
#  [5] "RNA_snn_res.0.8"  "seurat_clusters"  "CTgene"           "CTnt"
#  [9] "CTaa"   "CTstrict"    "clonalProportion" "clonalFrequency"
# [13] "cloneSize"
# Define color palette
colorblind_vector <- hcl.colors(
  n = 7, palette = "inferno", fixup = TRUE
)
colorblind_vector
# [1] "#040404" "#3E134F" "#851170" "#C53270" "#F36E35" "#F8B83C" "#FFFE9E"
DimPlot(
  seu_tcr,
  group.by = "cloneSize"
) +
  scale_color_manual(values = rev(colorblind_vector[c(1, 3, 5, 7)]))
