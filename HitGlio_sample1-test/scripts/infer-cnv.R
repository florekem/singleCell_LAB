# filtered counts matrix file
expr_matrix <- as.matrix(
  GetAssayData(
    seu_singlet,
    layer = "counts", assay = "SCT"
  )
)

expression_df <- as.data.frame(expr_matrix)
dim(expression_df)
head(expression_df)


write.table(expression_df, file = "expression_matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

test1 <- read_tsv("expression_matrix.txt")
dim(test1)

# sample annotation file
sample_annotation <- seu_singlet@meta.data |>
  rownames_to_column("cells") |>
  as_tibble() |>
  # select(cells, cluster.totalvi)
  select(cells, wsnn_res.0.4)

write.table(sample_annotation, file = "sample_annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

test2 <- read_tsv("sample_annotation.txt")
dim(test2)
test2
head(test1)

# gene ordering file
# "gene_ordering.txt"
# downloaded from: https://data.broadinstitute.org/Trinity/CTAT/cnv/

# run inferCNV
library(infercnv)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = "expression_matrix.txt",
  annotations_file = "sample_annotation.txt",
  delim = "\t",
  gene_order_file = "gene_ordering.txt",
  ref_group_names = NULL
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1, # may be 0 in case if filtering done in Seurat
  out_dir = "inferCNV_output_wsnn04", # dir is auto-created for storing outputs
  cluster_by_groups = T, # spr. jak wyjdzie na F
  denoise = T,
  HMM = T,
  # leiden_resolution = 0.01,
  num_threads = 40
)

?infercnv::run


str(expr_matrix)


dim(expr_matrix)
dim(sample_annotation)
