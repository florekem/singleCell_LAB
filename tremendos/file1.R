library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(tidyverse)

setwd("/mnt/sda4/singleCell_LAB/tremendos")

resOrdered <- readRDS("resOrdered.rds")
resOrdered

gsea_tumor_c3 <- readRDS("resOrdered_ESvsEF.rds")
length(gsea_tumor_c3)
upreg <- gsea_tumor_c3 > 0
upreg_gsea_tumor_c3 <- gsea_tumor_c3[upreg]
downreg <- gsea_tumor_c3 < 0
downreg_gsea_tumor_c3 <- gsea_tumor_c3[downreg]
length(upreg_gsea_tumor_c3)
length(downreg_gsea_tumor_c3)

ggo <- groupGO(
  gene = names(gsea_tumor_c3),
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  level = 3,
  readable = TRUE
)
head(ggo)

ego <- enrichGO(
  gene = names(downreg_gsea_tumor_c3),
  universe = resOrdered,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "MF",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)
dotplot(ego)

ego3_tumor_c3 <- gseGO(
  geneList = gsea_tumor_c3,
  OrgDb = "org.Mm.eg.db",
  ont = "BP",
  keyType = "SYMBOL",
  minGSSize = 2,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)
dotplot(ego3_tumor_c3, showCategory = 30, x = "enrichmentScore")



# over-representation analysis of ANY (not only go) database (msig)
# msig files downloaded from the website.
# i've tried to use msigdbr package, but it is outdated
# so now I can use msigdb not only for gsea, but also for over-rep.
msig <- read.gmt("mh.all.v2024.1.Mm.symbols.gmt")
msig <- read.gmt("./m2.all.v2024.1.Mm.symbols.gmt")
msig <- read.gmt("./m8.all.v2024.1.Mm.symbols.gmt")
colnames(msig) <- c("gs_name", "gene_symbol")

unique(msig$gs_name)


em <- enricher(
  names(upreg_gsea_tumor_c3),
  TERM2GENE = msig,
  universe = resOrdered,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
dotplot(em)

# there is not too much hits.
# now I want to extract genes for pathways, that contains specific
# keyword, and make heatmaps using theese genes, even if those pathways
# are not significant.
# search in msigdb .gmt for specific pathways and use them to
# construct heatmaps of genes that constitute those pathways
# creates a .rds file with list of names lists,
# where each name is a name of a pathway
msig_filtered <- filter(msig, grepl("", gs_name, ignore.case = TRUE))
msig_gs_name <- unique(msig_filtered$gs_name)
msig_gs_name <- as.character(msig_gs_name)
list_of_pathways_w_genes <- list()
for (i in msig_gs_name) {
  print(i)
  pathway <- msig_filtered |>
    filter(grepl(i, gs_name)) |>
    select(gene_symbol) |>
    as.vector()
  names(pathway) <- i
  list_of_pathways_w_genes <- append(list_of_pathways_w_genes, pathway)
}
list_of_pathways_w_genes
saveRDS(list_of_pathways_w_genes, "m2.hallmarl_all_pathways.rds")





em2 <- GSEA(
  gsea_tumor_c3,
  TERM2GENE = msig
)
