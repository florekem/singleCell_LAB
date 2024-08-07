# first integrate ADT assay, than RNA assay.
# The last integraion is used than in vizualization.
# Integrate ADT ------------
# Integrate ADT (może nie trzeba tego wcale robić?) ------------
# jako że joinowałem tylko assay hto, mogę bezkarnie integrować, bez dzielenia
DefaultAssay(seu_joined_singlet) <- "ADT"

seu_joined_singlet <- NormalizeData(
    seu_joined_singlet,
    normalization.method = "CLR"
)
seu_joined_singlet <- FindVariableFeatures(seu_joined_singlet, dims = 1:3)
seu_joined_singlet <- ScaleData(seu_joined_singlet)
seu_joined_singlet <- RunPCA(seu_joined_singlet)
seu_joined_singlet <- FindNeighbors(
    seu_joined_singlet,
    reduction = "pca",
    dims = 1:3
)
seu_joined_singlet <- FindClusters(
    seu_joined_singlet,
    resolution = 0.3, verbose = FALSE
)
seu_joined_singlet <- RunUMAP(
    seu_joined_singlet,
    reduction = "pca", dims = 1:3
)

seu_joined_singlet <- IntegrateLayers(
    object = seu_joined_singlet, method = CCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.adt.cca",
    assay = "ADT",
    verbose = FALSE
)
seu_joined_singlet <- FindNeighbors(
    seu_joined_singlet,
    reduction = "pca", dims = 1:3
)
seu_joined_singlet <- FindClusters(
    seu_joined_singlet,
    resolution = 0.3, verbose = FALSE
)
