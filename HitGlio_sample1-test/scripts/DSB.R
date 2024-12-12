library(dsb)
library(Seurat)
library(ggplot2)

setwd("/mnt/sda4/singleCell_LAB/HitGlio_sample1-test")

source("scripts/functions.R")

sample_path <- file.path("/mnt/sdb1/runs/sample1_multilane_spec_index_NNN") # nolint

# read raw data using the Seurat function "Read10X"
raw <- Seurat::Read10X(
  file.path(
    sample_path,
    "outs/multi/count/raw_feature_bc_matrix/"
  )
)
cells <- Seurat::Read10X(
  file.path(
    sample_path,
    "outs/per_sample_outs/sample1/count/sample_filtered_feature_bc_matrix/"
  )
)

# define cell-containing barcodes and separate cells and empty drops
stained_cells <- colnames(cells$`Gene Expression`)
background <- setdiff(colnames(raw$`Gene Expression`), stained_cells)

# split the data into separate matrices for RNA and ADT
prot <- raw$`Antibody Capture`
rna <- raw$`Gene Expression`

# create metadata of droplet QC stats used in standard scRNAseq processing
mtgene <- grep(pattern = "^MT-", rownames(rna), value = TRUE) # used below

md <- data.frame(
  rna.size = log10(Matrix::colSums(rna)),
  prot.size = log10(Matrix::colSums(prot)),
  n.gene = Matrix::colSums(rna > 0),
  mt.prop = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
)

# add indicator for barcodes Cell Ranger called as cells
md$drop.class <- ifelse(rownames(md) %in% stained_cells, "cell", "background")

# remove barcodes with no evidence of capture in the experiment
md <- md[md$rna.size > 0 & md$prot.size > 0, ]

ggplot(md, aes(x = log10(n.gene), y = prot.size)) +
  theme_bw() +
  geom_bin2d(bins = 300) +
  scale_fill_viridis_c(option = "C") +
  facet_wrap(~drop.class)

background_drops <- rownames(
  md[md$prot.size > 1.5 &
    md$prot.size < 3 &
    md$rna.size < 2.5, ]
)
background.adt.mtx <- as.matrix(prot[, background_drops])
dim(background.adt.mtx)

# calculate statistical thresholds for droplet filtering.
cellmd <- md[md$drop.class == "cell", ]

# filter drops with + / - 3 median absolute deviations from the median library size
rna.mult <- (3 * mad(cellmd$rna.size))
prot.mult <- (3 * mad(cellmd$prot.size))
rna.lower <- median(cellmd$rna.size) - rna.mult
rna.upper <- median(cellmd$rna.size) + rna.mult
prot.lower <- median(cellmd$prot.size) - prot.mult
prot.upper <- median(cellmd$prot.size) + prot.mult

# filter rows based on droplet qualty control metrics
qc_cells <- rownames(
  cellmd[cellmd$prot.size > prot.lower &
    cellmd$prot.size < prot.upper &
    cellmd$rna.size > rna.lower &
    cellmd$rna.size < rna.upper &
    cellmd$mt.prop < 0.14, ]
)

length(qc_cells)

# filter
cell.adt.raw <- as.matrix(prot[, qc_cells])
cell.rna.raw <- rna[, qc_cells]
cellmd <- cellmd[qc_cells, ]

pm <- sort(apply(cell.adt.raw, 1, max))
pm2 <- apply(background.adt.mtx, 1, max)
pm

# define isotype controls
isotype.controls <- c(
  "TotalSeqC_IgG2a_Ctrl",
  "TotalSeqC_IgG2b_Ctrl",
  "TotalSeqC_IgG1_Ctrl",
  "TotalSeqC_IgG_Ctrl"
)

# normalize and denoise with dsb with
cells.dsb.norm <- DSBNormalizeProtein(
  cell_protein_matrix = cell.adt.raw,
  empty_drop_matrix = background.adt.mtx,
  denoise.counts = TRUE,
  use.isotype.control = TRUE,
  isotype.control.name.vec = isotype.controls
)

# Seurat workflow
# integrating with Seurat
stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.adt.raw))))
stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.rna.raw))))

# create Seurat object note: min.cells is a gene filter, not a cell filter
s <- Seurat::CreateSeuratObject(
  counts = cell.rna.raw,
  meta.data = cellmd,
  assay = "RNA",
  min.cells = 20
)

# add dsb normalized matrix "cell.adt.dsb" to the "CITE" data (not counts!) slot
s[["CITE"]] <- Seurat::CreateAssayObject(data = cells.dsb.norm)

# define proteins to use in clustering (non-isptype controls)
prots <- rownames(s@assays$CITE@data)[3:40]

# cluster and run umap
s <- Seurat::FindNeighbors(
  object = s, dims = NULL, assay = "CITE",
  features = prots, k.param = 30,
  verbose = FALSE
)

# direct graph clustering
s <- Seurat::FindClusters(
  object = s, resolution = 1,
  algorithm = 3,
  graph.name = "CITE_snn",
  verbose = FALSE
)
# umap (optional)
# s = Seurat::RunUMAP(object = s, assay = "CITE", features = prots,
#                     seed.use = 1990, min.dist = 0.2, n.neighbors = 30,
#                     verbose = FALSE)

# make results dataframe
d <- cbind(
  s@meta.data,
  as.data.frame(t(s@assays$CITE@data))
  # s@reductions$umap@cell.embeddings)
)

library(magrittr)
# calculate the median protein expression separately for each cluster
adt_plot <- d %>%
  dplyr::group_by(CITE_snn_res.1) %>%
  dplyr::summarize_at(.vars = prots, .funs = median) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("CITE_snn_res.1")

# plot a heatmap of the average dsb normalized values for each cluster
pheatmap::pheatmap(t(adt_plot),
  color = viridis::viridis(25, option = "B"),
  fontsize_row = 8, border_color = NA
)
