# slingshot

library(slingshot)

# Extract PCA embeddings and cluster labels
pca <- Embeddings(seu_singlet_myelo, "pca.myelo.sct")
class(pca)
head(pca)

clusters <- seu_singlet_myelo$cluster.sct.myelo
class(clusters)

# Run Slingshot
slingshot_obj <- slingshot(pca, clusterLabels = clusters)

# Plot trajectory
plot(slingshot_obj, col = clusters)
plot(slingshot_obj, main = "Trajectory in PCA Space", type = "lineages")


summary(slingshot_obj)

slingshot_obj
head(slingshot_obj$slingPseudotime_1)


pseudotime <- slingPseudotime(slingshot_obj)
head(pseudotime)

seu_singlet_myelo$pseudotime <- pseudotime[, 1]

DimPlot(seu_singlet_myelo, reduction = "umap.myelo.pca.sct", group.by = "pseudotime")



# monocle3 + pyEvoCell
library(monocle3)

# Convert SingleCellExperiment object to CellDataSet
seu_singlet_myelo
cds <- SeuratWrappers::as.cell_data_set(seu_singlet_myelo)
cds
## Warning: Monocle 3 trajectories require cluster partitions, which Seurat does not calculate. Please run 'cluster_cells' on your cell_data_set object

# cds@colData$seurat_clusters <- seu_singlet_myelo$cluster.sct.myelo

# cds <- preprocess_cds(cds, num_dim = 10)

## Step 3: Reduce the dimensions using UMAP
# cds <- reduce_dimension(cds, reduction_method = "UMAP")

# plot_cells(
#   cds,
#   color_cells_by = "seurat_clusters",
#   show_trajectory_graph = FALSE
# )


# Extract UMAP embeddings from Seurat
# umap_coords <- Embeddings(seu_singlet_myelo, reduction = "umap.myelo.pca.sct")
umap_coords <- Embeddings(seu_singlet_myelo, reduction = "wnn.umap")

# Add UMAP embeddings to Monocle3
reducedDims(cds)$UMAP <- umap_coords

## Step 4: Cluster the cells
cds <- cluster_cells(cds, reduction_method = "UMAP")

# plot_cells(cds, show_trajectory_graph = FALSE)

cds <- learn_graph(cds, use_partition = TRUE)

plot_cells(cds,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  color_cells_by = "seurat_clusters"
)

cds <- order_cells(cds)

plot_cells(cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 1.5,
)

# save monocle object
save_monocle_objects(cds = cds, directory_path = "/mnt/sda4/monocle_temp/my_cds_objects")


# Progressions
progressions <- data.frame(
  cell_id = rownames(cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex),
  milestone_id = cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex[, 1]
)
write.csv(progressions, "progressions.csv", row.names = FALSE)

# Milestone percentages
milestone_percentages <- data.frame(
  cell_id = rownames(principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex),
  milestone_id = principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex[, 1],
  percentage = 1 # Placeholder percentage (can be replaced with real values)
)
write.csv(milestone_percentages, "milestone_percentages.csv", row.names = FALSE)

# Extract UMAP coordinates
dimred <- as.data.frame(reducedDims(cds)$UMAP)
dimred$cell_id <- rownames(dimred)
write.csv(dimred, "dimred.csv", row.names = FALSE)

# Extract trajectory edges
edges <- as.data.frame(as.matrix(principal_graph(cds)$UMAP))
write.csv(edges, "trajectory_edges.csv", row.names = FALSE)

# Extract metadata
metadata <- as.data.frame(colData(cds))
# Add cell_id as a column
metadata$cell_id <- rownames(metadata)
# Save metadata to CSV
write.csv(metadata, "metadata.csv", row.names = FALSE)

# Extract count matrix
count_data <- as.data.frame(as.matrix(assay(cds)))
write.csv(count_data, "count_data.csv", row.names = TRUE)




gene_fits <- fit_models(cds, model_formula_str = "~cluster")


gene_module_df <- find_gene_modules(cds, resolution = 1e-2)
gene_module_df

plot_cells(
  cds,
  genes = gene_module_df,
  show_trajectory_graph = FALSE
)

cds@reduce_dim_aux
