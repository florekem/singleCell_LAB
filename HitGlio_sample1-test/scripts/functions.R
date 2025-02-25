mycolor2 <- viridis::scale_color_viridis(option = "viridis")
# Create sparce matricies
seu_create_sparse_matrix_list <- function(sample_path) {
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
  length(joint_bcs)

  gex <- gex[, joint_bcs]
  adt <- adt[, joint_bcs]

  hto <- as.matrix(hto[, joint_bcs]) # not sure if now or later

  sparse_matricies_list <- list("gex" = gex, "adt" = adt, "hto" = hto)

  return(sparse_matricies_list)
}
# Identify empty cells with DropletUtils
filter_gex_matrix_dropletutils <- function(sparse_matricies_list) {
  gex <- sparse_matricies_list$gex
  hto <- sparse_matricies_list$hto
  adt <- sparse_matricies_list$adt

  e_out <- DropletUtils::emptyDrops(gex)
  is_cell <- e_out$FDR <= 0.01
  table(is_cell)
  is_cell[is.na(is_cell)] <- FALSE
  gex_filt <- gex[, is_cell]

  stained_cells <- colnames(gex_filt)
  hto_filt <- hto[, stained_cells]
  adt_filt <- adt[, stained_cells]

  filtered_sparse_matricies_list <- list(
    "gex_filt" = gex_filt, "adt_filt" = adt_filt, "hto_filt" = hto_filt
  )

  return(filtered_sparse_matricies_list)
}
# Create Seurat Object
create_seurat_object <- function(filtered_sparse_matricies_list) {
  gex_filt <- filtered_sparse_matricies_list$gex_filt
  adt_filt <- filtered_sparse_matricies_list$adt_filt
  hto_filt <- filtered_sparse_matricies_list$hto_filt

  seu <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(as.matrix(gex_filt), sparse = TRUE)
  )
  seu <- Seurat::NormalizeData(seu)
  seu <- Seurat::FindVariableFeatures(seu, selection.method = "mean.var.plot")
  seu <- Seurat::ScaleData(seu, features = Seurat::VariableFeatures(seu))

  seu[["HTO"]] <- Seurat::CreateAssayObject(counts = hto_filt)
  seu[["ADT"]] <- Seurat::CreateAssayObject(counts = adt_filt)

  seu <- Seurat::NormalizeData(seu, assay = "HTO", normalization.method = "CLR")
  seu <- Seurat::NormalizeData(seu, assay = "ADT", normalization.method = "CLR")

  return(seu)
}

# Demultiplex
seu_demultiplex <- function(seu) {
  seu_demuxed <- Seurat::MULTIseqDemux(
    seu,
    assay = "HTO",
    quantile = 0.7,
    autoThresh = TRUE,
    maxiter = 5,
    qrange = seq(from = 0.1, to = 0.9, by = 0.05),
    verbose = TRUE
  )
  table(seu_demuxed$MULTI_ID)

  return(seu_demuxed)
}
# Subset singlets
subset_singlets <- function(seu_demuxed) {
  seu_singlet <- subset(seu_demuxed, idents = c("tumor", "csf"))
  seu_negative <- subset(seu_demuxed, idents = "Negative")

  return(list("seu_singlet" = seu_singlet, "seu_negative" = seu_negative))
}
# Add quality features to singets
add_quality_features <- function(seu) {
  seu <- Seurat::PercentageFeatureSet(
    seu,
    pattern = "^RP[SL]",
    col.name = "percent_ribo"
  )
  seu <- Seurat::PercentageFeatureSet(
    seu,
    pattern = "^MT-",
    col.name = "percent_mito"
  )
  seu <- Seurat::PercentageFeatureSet(
    seu,
    pattern = "^HB[^(P)]",
    col.name = "percent_hemoglob"
  )
  seu[["log10GenesPerUmi"]] <-
    log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)

  return(seu)
}
# Add cell cycle scoring
add_cell_cycle_scoring <- function(seu_singlet) {
  s_genes <- Seurat::cc.genes$s.g
  g2m_genes <- Seurat::cc.genes$g2m.genes

  seu_singlet <- Seurat::CellCycleScoring(
    seu_singlet,
    s.features = s_genes,
    g2m.features = g2m_genes,
    set.ident = TRUE
  )

  return(seu_singlet)
}

remove_doublets_w_dfinder <- function(seu_singlet) {
  homotypic_prop <- DoubletFinder::modelHomotypic(seu_singlet$cluster.sct)
  nExp_poi <- round(0.075 * nrow(seu_singlet)) # estimate ~7.5% doublets
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))

  seu_singlet <- DoubletFinder::doubletFinder(
    seu_singlet,
    PCs = 1:10,
    pN = 0.25,
    pK = 0.09,
    nExp = nExp_poi_adj,
    reuse.pANN = FALSE,
    sct = TRUE
  )
  colnames(seu_singlet@meta.data)
  head(seu_singlet@meta.data)

  return(seu_singlet)
}

# function for filtering out genes not present in my suerat object,
# required for Nebulosa::plot_density() as it required all genes to be
# present in the dateset
# could be replaced in-place by map(), for example:
# function(x) x[x %in% rownames(seu_singlet_myelo)]
valid_features <- function(features_list, seurat_obj) {
  features_list[features_list %in% rownames(seurat_obj)]
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

seu_gsva <- function(seurat_object, category, gsva_type) {
  gene_sets <- msigdbr::msigdbr(
    species = "Homo sapiens",
    category = category
  )
  gene_set_list <- split(gene_sets$gene_symbol, gene_sets$gs_name)
  expr_matrix <- as.matrix(
    Seurat::GetAssayData(
      seurat_object,
      layer = "data", assay = "RNA"
    )
  )

  if (gsva_type == "gsva") {
    params <- GSVA::gsvaParam(
      expr_matrix,
      gene_set_list,
      minSize = 5,
      maxSize = 500
    )
  } else if (gsva_type == "zscore") {
    params <- GSVA::zscoreParam(
      expr_matrix,
      gene_set_list,
      minSize = 5,
      maxSize = 500
    )
  } else if (gsva_type == "plage") {
    params <- GSVA::plageParam(
      expr_matrix,
      gene_set_list,
      minSize = 5,
      maxSize = 500
    )
  }

  gsva_es <- GSVA::gsva(
    params,
    verbose = FALSE
  )

  return(gsva_es)
}

seu_gsva_print_plot <- function(top_markers, title) {
  for (i in unique(top_markers$cluster)) {
    top_cluster <- top_markers |>
      filter(cluster == i) |>
      arrange(desc(avg_log2FC))

    print(
      ggplot2::ggplot(
        top_cluster, aes(x = gene, cluster, y = avg_log2FC, fill = upordown)
      ) +
        geom_bar(stat = "identity", width = 0.4) +
        coord_flip() +
        theme_minimal() +
        scale_fill_manual(values = c("#555598", "#bf605b")) +
        theme(
          axis.text = element_text(size = 11)
        ) +
        patchwork::plot_annotation(
          title = paste0(title, i),
          theme = theme(plot.title = element_text(size = 25))
        )
    )
  }
}

seu_gsea <- function(ranked_genes, category) {
  gene_sets <- msigdbr::msigdbr(
    species = "Homo sapiens",
    category = category
  )
  gene_set_list <- split(gene_sets$gene_symbol, gene_sets$gs_name)

  fgsea_results <- fgsea(
    pathways = gene_set_list,
    stats = ranked_genes,
    minSize = 5, maxSize = 500
  )

  top_pathways <- fgsea_results %>%
    filter(padj < 0.05) |>
    arrange(desc(NES)) |>
    slice_head(n = 20) |> # Top 20 positive
    bind_rows(
      fgsea_results |>
        filter(padj < 0.05) |>
        arrange(NES) |>
        slice_head(n = 20)
    )

  return(top_pathways)
}

seu_gsea_plot <- function(top_pathways, category, cluster) {
  plot1 <- ggplot(top_pathways, aes(reorder(pathway, NES), NES, fill = padj < 0.05)) +
    geom_col() +
    coord_flip() +
    labs(x = "Pathway", y = "NES", title = paste0("GSEA Top 40(20up+20down) ", category, " pathways in cluster ", cluster))

  print(plot1)
}



# redundant, not usign it, keeping for the future
seu_ssgsea <- function(seurat_object, collections_list, top_n) {
  GS.hallmark <- escape::getGeneSets(library = collection)

  seurat_object <- escape::runEscape(
    seurat_object,
    method = "ssGSEA",
    gene.sets = GS.hallmark,
    groups = 5000,
    min.size = 5,
    new.assay.name = "escape.ssGSEA"
  )

  plot <- escape::heatmapEnrichment(
    seurat_object,
    group.by = "manual_anno",
    gene.set.use = rownames(seurat_object@assays$escape.ssGSEA)[1:top_n],
    assay = "escape.ssGSEA",
    scale = TRUE,
    cluster.rows = TRUE,
    cluster.columns = TRUE
  ) +
    ggplot2::theme(
      # axis.text = element_text(size = 11),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
    )

  gsea_plots <- c(gsea_plots, list(plot))

  return(gsea_plots)
}
