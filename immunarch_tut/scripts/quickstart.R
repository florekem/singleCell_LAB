library(immunarch)
data(immdata)

class(immdata)

repExplore(immdata$data, "lens") %>% vis()


data(scdata)
