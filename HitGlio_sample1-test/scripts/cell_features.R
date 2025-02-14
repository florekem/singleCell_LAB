# NOLINTBEGIN

# /random
f1 <- c("CD8B", "NKG7")
f2 <- c("CD4")
f3 <- c("CX3CR1")
f4 <- c("H3F3A", "IDH1", "IDH2", "SETD2")
f5 <- c("FOXP3", "TOX", "CD4")

# /cell metaprograms /
# /from: Liu, Ilon, Li Jiang, Erik R. Samuelsson, Sergio Marco Salas, Alexander Beck, Olivia A. Hack, Daeun Jeong, et al. “The Landscape of Tumor Cell States and Spatial Organization in H3-K27M Mutant Diffuse Midline Glioma across Age and Location.” Nature Genetics 54, no. 12 (December 2022): 1881–94. https://doi.org/10.1038/s41588-022-01236-3. /# nolint

liu2022 <- list(
  ac_like = c("CLU", "AQP4", "AGT", "SPARCL1", "VIM", "GFAP", "CRYAB", "MLC1", "PLTP", "LAMB2", "CD99", "APOE", "S1PR1", "SPARC", "GJA1", "EDNRB", "SLC1A3", "F3", "CD44", "TNC", "C4A", "EZR", "HEPN1", "ID3", "CCDC80", "S100A10", "IFI6", "RASSF4", "ATP1A2", "ALDOC"), # nolint
  mes_like = c("VIM", "S100A10", "TIMP1", "GAP43", "CLU", "S100A6", "TAGLN2", "TMSB10", "MYL12A", "S100A16", "DDIT3", "SPP1", "LMNA", "TNFRSF12A", "VMP1", "IL13RA2", "SQSTM1", "SLC2A3", "OCIAD2", "CD63", "SPARCL1", "LGALS1", "MT2A", "SCG2", "S100A11", "CAV1", "CD59", "EMP1"),
  opc_like1 = c(
    "PDGFRA", "APOD", "PLLP", "EPN2", "MT3", "COL9A1", "CA10", "RAMP1", "FABP7", "OLIG1",
    "MEG3", "LPPR1", "OMG", "CHAD", "TSPAN7", "HSPA1B", "LHFPL3", "CD9", "ACSF2", "GPR17",
    "NSG2", "S100B", "KLRC2", "RP11-231C18.1", "CST3", "TM4SF1", "ABHD2", "TRIO", "NKAIN4",
    "RP11-87N24.3"
  ),
  opc_like2 = c(
    "RPS4Y1", "CRABP1", "RPS17", "RPS21", "RPS18", "EPB41L4A-AS1", "RPL12", "RPL34", "PCP4",
    "RPS25", "RPL22L1", "SNHG8", "RPL36", "RPS14", "ZFAS1", "RPS17L", "RPL26", "RPL36A",
    "RPS29", "RPL27", "RPL31", "RPL23", "PRAME", "RPS27", "RPL38", "NR4A1", "BTG2",
    "RPL17", "RPS12", "RPS13"
  ),
  opc_like3 = c("IER3", "CTC-250I14.6", "FOS", "KLF6", "JUNB", "ZFP36", "BTG2", "ARC", "ATF3", "NFKBIZ", "DLL1", "DUSP10", "RFX4", "RP11-698N11.2", "FOSB", "GADD45B", "IFI16", "XIST", "EGR2", "CYR61", "DUSP1", "MEST", "EGFR", "EGR1", "LINC00910", "ID2", "JUN", "IRF1", "PRCP", "CREB5"),
  oc_like = c("PTGDS", "RP11-862L9.3", "MBP", "PLP1", "TF", "BCAS1", "MOG", "SIRT2", "TUBB4A", "GPR17", "CLDN11", "APLP1", "CNP", "CDK18", "MYRF", "MAG", "RGR", "UGT8", "MPZL1", "NFASC", "LPPR1", "FYN", "GNAI1", "ZNF488", "BMPER", "DBNDD2", "FRMD4B", "SEMA4D", "FAM13C", "SGK1"),
  s = c("KIAA0101", "PCNA", "TYMS", "RRM2", "MLF1IP", "UBE2T", "MCM4", "ZWINT", "CDK1", "MAD2L1", "GMNN", "HMGB2", "SMC4", "DHFR", "CENPK", "MCM5", "FANCI", "DTL", "HELLS", "EXO1", "DUT", "MCM2", "TK1", "RAD51AP1", "CENPM", "PBK", "TOP2A", "LIG1", "FEN1", "PKMYT1"),
  g2m = c("HMGB2", "UBE2C", "TOP2A", "CDK1", "UBE2T", "NUSAP1", "PRC1", "MAD2L1", "PTTG1", "PBK", "CCNB1", "TYMS", "FAM64A", "SMC4", "RRM2", "CENPF", "NUF2", "CCNB2", "ZWINT", "TPX2", "RACGAP1", "CKS2", "MLF1IP", "BIRC5", "CDC20", "BUB1B", "CDKN3", "SPAG5", "KPNA2", "KIAA0101")
)
# nolint
# /mes like not available in out dataset/
# /"RP11-430H10.1"/
# /"ORAOV1"/

# /cell types from: Neftel, Cyril, Julie Laffy, Mariella G. Filbin, Toshiro Hara, Marni E. Shore, Gilbert J. Rahme, Alyssa R. Richman, et al. “An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma.” Cell 178, no. 4 (August 8, 2019): 835-849.e21. https://doi.org/10.1016/j.cell.2019.06.024. / #nolint
oligoden <- c("MBP", "TF", "PLP1", "MAG", "MOG", "CLDN11")
tcels <- c("CD2", "CD3D", "CD3E", "CD3G")
macrophages <- c("CD14", "AIF1", "FCER1G", "FCGR3A", "TYROBP", "CSF1R")
t_cell_markers <- c(
  "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B",
  "TRAC", "TRBC1", "TRBC2", "IL7R", "FOXP3",
  "CD25", "CD69", "CD28", "CTLA4", "PDCD1",
  "CCR7", "SELL", "GZMB", "PRF1", "IFNG"
)

ein_like <- c(
  "MEG3", "STMN2", "DLX6-AS1", "RTN1", "TUBB3", "GAD1", "NSG1", "UCHL1", "PLS3",
  "RUNX1T1", "DLX5", "NRXN3", "MAP1B", "TMSB10", "ARL4D", "STMN4", "CADM1",
  "ERBB4", "DLX6", "FAM65B", "MIAT", "LMO3", "SLAIN1", "TPGS2", "BCL11A",
  "CACNA2D3", "NCAM1", "GNG2", "AC017053.1", "CRMP1", "GAP43", "ATCAY",
  "TP53I11", "PLK2", "WLS", "RP11-87N24.3", "RNF24", "STXBP1", "MYCBP2",
  "EPHA3", "MEIS2", "TMSB15A", "ENO2", "LMBR1L", "NSG2", "GPD1", "SCAPER",
  "CELF4", "KLF7", "TMEM123"
)
inpc_like <- c(
  "MFNG", "TAC3", "INPPL1", "BEST3", "DLL1", "CALCOCO1", "CBFA2T2", "HDAC9", "RGS8",
  "RP11-87N24.3", "DLL3", "RP11-118F19.1", "FTX", "JAG1", "TCF12", "ENC1",
  "ARHGEF2", "TMEM2", "MAP6", "ITFG1", "LAMA5", "COPA", "ANGPTL2", "AKNAD1",
  "LAMP5", "DLGAP1", "DLG2", "BOC", "KAT7", "KIAA0195", "EIF2A", "HIST1H2AC",
  "NPY", "GLYCTK", "FAM60A", "SS18", "RGS16", "GPR56", "LIMD1", "CNTN2",
  "JAM3", "CORO7", "QARS", "AC013402.2", "CDK5RAP1", "ARHGAP21", "ACADVL",
  "THOC1", "ADAM10", "PXDN"
)


tam_gene_sets <- list( # deepseek R1
  # Cell Types
  Monocytes = c("CD14", "FCGR3A", "LYZ", "S100A8", "S100A9"),
  M_MDSC = c("CD14", "S100A8", "S100A9", "ARG1", "IL4R"),
  E_MDSC = c("FUT4", "CEACAM8", "MPO", "STAT3"),
  Macrophages = c("CD68", "CD163", "CSF1R", "MRC1"),
  M1_TAM = c("IL1B", "TNF", "CD80", "CXCL10"),
  M2_TAM = c("MRC1", "CD163", "CCL18", "VEGFA"),
  LAM_TAM = c("TREM2", "APOE", "LPL", "LIPA"),

  # Functional States
  Hypoxia = c("HIF1A", "VEGFA", "CA9", "SLC2A1"),
  Chemotaxis = c("CXCL8", "CCL2", "CXCR4", "CCR2"),
  Angiogenesis = c("VEGFA", "ANGPT1", "PDGFB", "FLT1"),
  Immune_Suppression = c("CD274", "CTLA4", "IDO1", "IL10")
)



wang2024_hypoxic_macrophages <- list(
  # Cell Types
  Tumor_cell = c("GFAP", "PTPRZ1", "SCG3", "Clorf61", "SOX2", "SOX6", "BCAN", "NOVA1", "NRCAM", "TUBB2B"),
  Proliferating_tumor_cell = c("SOX2", "SOX6", "MK167", "TOP2A", "CCDC34", "PBK", "UBE2T", "CENPV", "CENPF"),
  Oligodendrocyte = c("MAG", "CLDN11", "APLP1", "TMEM144", "CNDP1", "EDIL3", "LARP6", "MOG", "AMER2", "TUBB4A"),
  NK_T = c("CD3D", "CD3E", "TRBC2", "CD3G", "CD2", "CD52", "TRAC", "ILTR", "TRBC1", "LCK", "SKAP1", "CD48"),
  Proliferating_NK_T = c("MK167", "TOP2A", "CD3D", "CD3E", "TRBC2", "CD3G"),
  Pericyte = c("COLL42", "COL6A2", "PDGFRB", "COL3A1", "EDNRA", "LUM", "COL1A1", "PLAC9", "FRZB"),
  Endothelial_cell = c("ESAM", "FLT1", "RAMP2", "VWF", "EGFL7", "SLC943R2", "CAVIN2", "ADGRL4", "ABLIM1", "ABCB1"),
  Dendritic_cell = c("AREG", "FCER1A", "CLEC10A", "HLA-DQB1", "HLA-DQA1"),
  Monocyte = c("FCN1", "VCAN", "CD52", "S100A8", "S100A9"),

  # Tumor-Associated Macrophage (TAM) States
  Proliferating_TAM = c("MK167", "STMN1", "TYMS", "TOP2A", "TUBB"),
  Microglia_TAM = c("P2RY12", "CX3CR1", "FCGR1A", "CH25H", "CCL4L2"),
  Monocyte_derived_TAM = c("TGFBI", "CD14", "CD163", "SELENOP", "GPNMB", "CSTB", "ADM", "BNIP3", "BNIP3L", "ENO2", "FAM162A", "RALA", "CYTIP", "ADAM8", "C15orf48", "CD109"),
  Hypoxia_TAM = c("ENO1", "SCD", "PLP2", "RAB42", "S100A10", "S100A6", "HK2", "SLC2A1", "CXCL8", "LG4LS1", "TIMP1", "PLIN2", "CTSL", "LDHA", "NDRG1", "HILPDA", "ERO1A", "NUPR1", "MT2A"),
  Chemotaxis_TAM = c("SPP1", "RGS16", "FSCN1", "HAMP", "BIN1", "IBSP"),
  IFN_TAM = c("IFIT1", "IF144L", "IFIT3", "IF16", "IF144", "FGL2", "HERC5", "IFIT2", "ISG15", "CXCL10", "MX1", "MX2"),
  Lipid_TAM = c("APOC1", "APOE", "CD63", "LG4LS3", "ACP5", "GCHFR", "OTOA", "PLA2G7", "SDS"),
  Phago_AP_TAM = c("CD83", "CD74", "HLA-DRA", "HLA-DQA1", "DDIT3", "THBS1", "CD69", "FOSB", "TNFAIP3", "ZFP36", "BAG3"),
  Ribosome_TAM = c("RPL13", "RPL19", "RPS15", "RPL28", "RPLP1", "PRL35")
)



greenwald2024 <- list(
  inflammatory_macrophages = c(
    "CCL2", "CCL3", "CCL4", "CCL4L2", "DNAJB1", "EGR1", "IL1B", "NR4A1", "CD83", "FOS",
    "BTG2", "CCL3L1", "CH25H", "CXCL10", "CYR61", "DUSP1", "GBP1", "HSPA1A", "HSPA1B",
    "HSPA6", "HSPH1", "IER2", "IER5", "IRF1", "JUN", "JUNB", "RGS1", "ZFP36", "AHRR",
    "EGR3", "CLPTM1L", "TAP1", "STAT1", "SERPINE1", "EXOC3-AS1", "GADD45B", "SEZ6L",
    "IGFBP3", "ATF3", "SLPI", "ROPN1L", "GPX3", "AZGP1", "WARS", "CCDC127", "G0S2",
    "PLAUR", "VGF", "CD14", "APLNR"
  ),
  non_malignant_macrophages = c(
    "C1QB", "C1QA", "FCER1G", "C1QC", "LAPTM5", "CD14", "CD74", "SRGN", "AIF1", "FCGR3A",
    "TYROBP", "HLA-DRB1", "CSF1R", "HLA-DRA", "TREM2", "APOC1", "C3", "HAMP", "ITGB2",
    "VSIG4", "ALOX5AP", "FCGBP", "SPP1", "CD163", "HLA-DPA1", "CD68", "S100A8", "HLA-DRB5",
    "MS4A6A", "RGS1", "RNASET2", "S100A9", "HLA-DPB1", "CTSS", "SLCO2B1", "CD53", "CTSC",
    "HLA-DMB", "CHI3L1", "HLA-DQB1", "S100A11", "SLC2A5", "A2M", "EFEMP1", "GPNMB",
    "HLA-DQA1", "LILRB4", "STAB1", "VAMP8", "IBSP"
  )
)


hara2021 <- list(
  cc_genes = c("CENPF", "UBE2C", "NUSAP1", "PTTG1", "HMGB2", "UBE2T", "BIRC5", "CDKN3", "PRC1", "MKI67", "CENPW", "AURKB", "TOP2A", "NUF2", "CCNB1", "MAD2L1", "GTSE1", "UBE2S", "CKS2", "CCNB2", "SMC4", "PBK", "TPX2", "CDCA3", "CDC20"),
  ac_genes = c("AGT", "TTYH1", "SPARCL1", "ATP1B2", "CST3", "SCRG1", "KANK1", "CRYAB", "CDC42EP4", "HEY1", "SAT1", "APOE", "NUPR1", "BAALC", "GPC1", "SEMA6D", "APOC1", "CA12", "TFPT", "NCALD", "IFIT2", "PJA2", "C1orf54", "ADGRG1", "DKK3"),
  mes_genes = c("CD44", "ANXA1", "LGALS3", "IFITM3", "BST2", "S100A10", "CDKN1A", "RAMP1", "TIMP4", "ATOX1", "GDF15", "CD9", "COL4A2", "MT1X", "SPARC", "S100A13", "S100A6", "IGFBP5", "GOLIM4", "PHLDA2", "LAMA4", "CD163L1", "MUC12", "CD164", "IFI27L2"),
  npc_genes = c("STMN2", "DCX", "DLX5", "SOX4", "GADD45G", "RND3", "SCG3", "SOX11", "MLLT11", "GDAP1L1", "ELAVL4", "TAGLN3", "NREP", "NFIA", "TCF4", "CD24", "PCP4", "STMN4", "MAP1B", "NFIB", "SNTG1", "BEX1", "BCL7A", "MAP2", "KIF5C")
)



# NOLINTEND
