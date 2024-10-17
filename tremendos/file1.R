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


library(msigdbr)
msigdbr_species()
?msigdbr

mouse_df <- msigdbr(species = "mouse")
unique(mouse_df$gs_cat)
# [1] "C3" "C2" "C8" "C6" "C7" "C4" "C5" "H"  "C1"
head(mouse_df, 3) %>% as.data.frame()

m_t2g <- msigdbr(species = "mouse", category = "C6") %>%
  dplyr::select(gs_name, gene_symbol)
head(m_t2g)

# over-representation analysis
em <- enricher(names(gsea_tumor_c3), TERM2GENE = m_t2g)
head(em)

msig <- read.gmt("mh.all.v2024.1.Mm.symbols.gmt")
msig <- read.gmt("m2.all.v2024.1.Mm.symbols.gmt")
colnames(msig) <- c("gs_name", "gene_symbol")

unique(msig$gs_name)

# search in msigdb .gmt for specific pathways and use them to
# construct heatmaps of genes that constitute those pathways
msig_filtered <- filter(msig, grepl("phagocytosis", gs_name, ignore.case = TRUE))
msig_filtered

pathay_wp_microglia_phagocyt <- msig_filtered |>
  filter(grepl("WP_MICROGLIA_PATHOGEN_PHAGOCYTOSIS_PATHWAY", gs_name)) |>
  select(gene_symbol) |>
  as.vector()

pathay_ractome_role_of_phoph_in_phago <- msig_filtered |>
  filter(grepl("REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS", gs_name)) |>
  select(gene_symbol) |>
  as.vector()

em <- enricher(
  names(upreg_gsea_tumor_c3),
  TERM2GENE = msig,
  universe = resOrdered,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
dotplot(em)

em2 <- GSEA(
  gsea_tumor_c3,
  TERM2GENE = msig
)
