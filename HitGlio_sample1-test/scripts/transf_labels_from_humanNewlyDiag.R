# / https://brainimmuneatlas.org/data_files/toDownload/filtered_feature_bc_matrix_HumanNewlyDiagnGBM.zip #nolint

# seu_singlet from 01_demux_filter.R

# 1. Paths ------------------------------------------------------
sample_path <- file.path("/mnt/sda4/brain_immune_atlas/humanNewlyDiagnosedGBM")
print(sample_path)

# 2. Create seurat object ---------------------------------------
sparse_matrix <- Seurat::Read10X(data.dir = sample_path)

seu <- CreateSeuratObject(
  counts = Matrix::Matrix(as.matrix(sparse_matrix), sparse = TRUE)
)
seu

annot <- read.csv("/mnt/sda4/brain_immune_atlas/annot_Human_ND_GBM_Full.csv")
colnames(annot)

## 2.1. filtering by annotation file --------
seu <- subset(
  seu,
  cells = annot$cell
)
seu

## 2.2. add metadata from anot file ---------------------
seu <- AddMetaData(
  seu,
  metadata = annot
)
colnames(seu@meta.data)
Idents(seu) <- "cluster"
unique(Idents(seu))

# 3.1 process new dataset (prepare to integrate) ------------------
# Integration with log-normalized data
# it does not work with SCT data?
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
ElbowPlot(seu, ndims = 50)

seu <- FindNeighbors(
  seu,
  reduction = "pca", dims = 1:50 # increased bc of SCT
)
seu <- FindClusters(
  seu,
  resolution = 0.2, verbose = FALSE,
  algorithm = 4 # leiden
)
seu <- RunUMAP(
  seu,
  reduction = "pca", dims = 1:50
)
seu
DimPlot(
  seu,
  group.by = "cluster", reduction = "umap",
  label = TRUE
)
# 4. Merge both datasets --------------------------
DefaultAssay(seu) <- "RNA"
DefaultAssay(seu_singlet) <- "RNA"

merged <- merge(seu_singlet, seu)
merged
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)

# 5. Integrate ----------------------
integrated <- IntegrateLayers(
  merged,
  method = HarmonyIntegration,
  # method = CCAIntegration,
  orig.reduction = "pca",
  # assay = "RNA",
  new.reduction = "integrated.harmony"
)
integrated <- FindNeighbors(
  integrated,
  reduction = "integrated.harmony",
  dims = 1:30
)
integrated <- FindClusters(integrated) # possible to use Leiden?

# 6. Transfer labels -------------------
transf.anchors <- FindTransferAnchors(
  reference = integrated,
  query = seu_singlet,
  dims = 1:30,
  reference.reduction = "pca"
)

predictions <- TransferData(
  anchorset = transf.anchors,
  refdata = integrated$cluster,
  dims = 1:30
)

seu_singlet <- AddMetaData(
  seu_singlet,
  metadata = predictions
)

DefaultAssay(seu_singlet) <- "RNA"
colnames(seu_singlet@meta.data)
DimPlot(
  seu_singlet,
  group.by = "predicted.id", reduction = "umap",
  label = TRUE
)
table(seu_singlet@meta.data$predicted.id)

?DimPlot
seu_singlet
