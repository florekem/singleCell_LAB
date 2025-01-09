# Visualize features
inferno <- viridis::scale_color_viridis(option = "inferno")

q_features <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "log10GenesPerUmi",
  "percent_mito",
  "percent_ribo",
  "percent_hemoglob"
)

gr_by_features <- c("MULTI_ID", "cluster.log")
gr_by_features_sct <- c("MULTI_ID", "cluster.sct")
gr_by_features_sct_myelo <- c("MULTI_ID", "cluster.sct.myelo")
gr_by_features_totalvi <- c("MULTI_ID", "cluster.totalvi")
t_quality <- c("CD8A", "CD8B", "CD4")
t_adt_quality <- c("TotalSeqC-CD8", "TotalSeqC-CD4")
isotype_controls <- c(
  "TotalSeqC-IgG2b-Ctrl",
  "TotalSeqC-IgG2a-Ctrl",
  "TotalSeqC-IgG-Ctrl",
  "TotalSeqC-IgG1-Ctrl"
)
