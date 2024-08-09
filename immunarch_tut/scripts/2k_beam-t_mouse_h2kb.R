library(immunarch)

file_path <- "/mnt/sda4/singleCell_LAB/immunarch_tut/data"
immdata_10x <- repLoad(file_path)

names(immdata_10x)
# [1] "data" "meta"
immdata_10x$meta
repExplore(immdata_10x$data, "lens") %>% vis()

exp_vol <- repExplore(immdata_10x$data, .method = "volume")
p1 <- vis(exp_vol)
p1
