# /random
f1 <- c("CD8B", "NKG7")
f2 <- c("CD4")
f3 <- c("CX3CR1")
f4 <- c("H3F3A", "IDH1", "IDH2", "SETD2")
f5 <- c("FOXP3", "TOX", "CD4")

# /cell metaprograms /
# /from: Liu, Ilon, Li Jiang, Erik R. Samuelsson, Sergio Marco Salas, Alexander Beck, Olivia A. Hack, Daeun Jeong, et al. “The Landscape of Tumor Cell States and Spatial Organization in H3-K27M Mutant Diffuse Midline Glioma across Age and Location.” Nature Genetics 54, no. 12 (December 2022): 1881–94. https://doi.org/10.1038/s41588-022-01236-3. /# nolint
ac_like <- c("CLU", "AQP4", "AGT", "SPARCL1", "VIM", "GFAP", "CRYAB", "MLC1", "PLTP", "LAMB2", "CD99", "APOE", "S1PR1", "SPARC", "GJA1", "EDNRB", "SLC1A3", "F3", "CD44", "TNC", "C4A", "EZR", "HEPN1", "ID3", "CCDC80", "S100A10", "IFI6", "RASSF4", "ATP1A2", "ALDOC") # nolint
mes_like <- c("VIM", "S100A10", "TIMP1", "GAP43", "CLU", "S100A6", "TAGLN2", "TMSB10", "MYL12A", "S100A16", "DDIT3", "SPP1", "LMNA", "TNFRSF12A", "VMP1", "IL13RA2", "SQSTM1", "SLC2A3", "OCIAD2", "CD63", "SPARCL1", "LGALS1", "MT2A", "SCG2", "S100A11", "CAV1", "CD59", "EMP1") # nolint
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
