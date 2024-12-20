mycolor2 <- viridis::scale_color_viridis(option = "viridis")
?scale_color_viridis

initial_steps <- function(sample_path) {
  sparse_matrix <- Seurat::Read10X(data.dir = sample_path)
  gex <- sparse_matrix$`Gene Expression`
  adt <- sparse_matrix$`Antibody Capture`[grepl(
    "TotalSeqC", rownames(sparse_matrix$`Antibody Capture`)
  ), ]
  hto <- sparse_matrix$`Antibody Capture`[grepl(
    "Hashtag", rownames(sparse_matrix$`Antibody Capture`)
  ), ]
  rownames(hto)
  rownames(hto) <- c("tumor", "csf")
  joint_bcs <- intersect(colnames(gex), colnames(hto))
  head(joint_bcs)

  gex <- gex[, joint_bcs]
  adt <- adt[, joint_bcs]

  hto <- as.matrix(hto[, joint_bcs])

  # adt <- adt[, joint_bcs]
  # ncol(adt)
  # hto <- hto[, joint_bcs]
  # ncol(hto)

  seu <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(as.matrix(gex), sparse = TRUE)
  )

  # Demux and subset singlets -----------------------------------
  seu <- Seurat::NormalizeData(seu)
  seu <- Seurat::FindVariableFeatures(seu, selection.method = "mean.var.plot")
  seu <- Seurat::ScaleData(seu, features = Seurat::VariableFeatures(seu))

  # Add HTO data as a new assay independent from RNA ----
  seu[["HTO"]] <- Seurat::CreateAssayObject(counts = hto)
  seu[["ADT"]] <- Seurat::CreateAssayObject(counts = adt)

  seu <- Seurat::NormalizeData(seu, assay = "HTO", normalization.method = "CLR")

  return(seu)
}

initial_demux <- function(seu) {
  seu <- Seurat::MULTIseqDemux(
    seu,
    assay = "HTO",
    quantile = 0.7,
    autoThresh = TRUE,
    maxiter = 5,
    qrange = seq(from = 0.1, to = 0.9, by = 0.05),
    verbose = TRUE
  )
  return(seu)
}

seu_subset <- function(seu, idents = NULL, subset = NULL) {
  subset(seu,
    idents = idents,
    subset = subset
  )
}

seu_violin_plot <- function(seu, features) {
  Seurat::VlnPlot(
    seu,
    features = features
  )
}

viz_quality <- function(
    seu, assay, q_features = NULL,
    reduction, gr_by_features = NULL) {
  #
  purrr::walk(gr_by_features, ~ seu_DimPlot(
    seu,
    reduction = reduction,
    group.by = .x
  ) |> print())

  purrr::walk(q_features, ~ seu_FeaturePlot(
    seu,
    assay = assay,
    features = .x,
    reduction = reduction,
    color = inferno,
    label = TRUE,
    max.cutoff = NA
  ) |> print())

  purrr::walk(q_features, ~ seu_violin_plot(seu, features = .x) |> print())
}

plot_scatter <- function(data) {
  Seurat::FeatureScatter(
    data,
    "nFeature_RNA",
    "nCount_ADT"
  )
}

seu_normalize_var_scale_pca <- function(
    seu, normalize, run_pca = TRUE, dims, pca.reduction.name, k.param = 20,
    cluster.name, umap.reduction.name, algorithm, resolution,
    group.singletons = TRUE) {
  Seurat::DefaultAssay(seu) <- "RNA"
  if (normalize == "log") {
    print("running lognorm...")
    seu <- Seurat::NormalizeData(seu)
    seu <- Seurat::FindVariableFeatures(seu)
    seu <- Seurat::ScaleData(seu)
  }
  if (normalize == "sct") {
    print("running SCT...")
    seu <- Seurat::SCTransform(seu, verbose = FALSE)
  }
  if (normalize == FALSE) {
    print("Running WITHOUT Normalization")
  }
  if (run_pca == TRUE) {
    seu <- Seurat::RunPCA(
      seu,
      reduction.name = pca.reduction.name
    )
  } else {
    print("NOT RUNNING PCA!") # in case other dim. reductions such as scVI
  }
  seu <- Seurat::FindNeighbors(
    seu,
    reduction = pca.reduction.name,
    dims = dims,
    k.param = k.param
  )
  seu <- Seurat::FindClusters(
    seu,
    cluster.name = cluster.name,
    algorithm = algorithm, # leiden
    resolution = resolution,
    group.singletons = group.singletons,
    verbose = FALSE
  )
  seu <- Seurat::RunUMAP(
    seu,
    reduction = pca.reduction.name,
    reduction.name = umap.reduction.name, # new name
    dims = dims
  )
}
?FindClusters
?RunUMAP

seu_DimPlot <- function(
    seu, group.by, reduction, label = TRUE, label.size = 4,
    show = TRUE, save_path = FALSE, ggwidth = NA, ggheight = NA) {
  plot <- Seurat::DimPlot(
    seu,
    group.by = group.by,
    reduction = reduction,
    label = label,
    repel = TRUE,
    label.size = label.size
  )
  if (save_path != FALSE) {
    ggplot2::ggsave(paste0(save_path, ".png"), plot, width = ggwidth, height = ggheight)
  }
  if (show == TRUE) {
    return(plot)
  }
}
seu_FeaturePlot <- function(
    seu, assay, features, reduction, label = TRUE,
    max.cutoff = NA, color = mycolor2, save_path = FALSE, show = TRUE,
    ggwidth = NA, ggheight = NA, ...) {
  Seurat::DefaultAssay(seu) <- assay

  plot <- Seurat::FeaturePlot(
    seu,
    features = features,
    reduction = reduction,
    label = label,
    max.cutoff = max.cutoff
  ) & color
  if (save_path != FALSE) {
    ggplot2::ggsave(paste0(save_path, ".png"), plot, width = ggwidth, height = ggheight)
  }
  if (show == TRUE) {
    return(plot)
  }
}
seu_SCT <- function(seu, verbose = FALSE) {
  Seurat::SCTransform(
    seu,
    verbose = FALSE
  )
}
?FeaturePlot
seu_dittoDotPlot <- function(seu, vars, group.by) {
  dittoSeq::dittoDotPlot(
    seu,
    vars = vars,
    group.by = group.by
  )
}

get_features <- function(seu, assay) {
  DefaultAssay(seu) <- assay
  features <- Features(seu)

  return(features)
}
