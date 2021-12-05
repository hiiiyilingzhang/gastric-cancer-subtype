# Downloaded from Source Data Fig.1

# 1. Read files ----
gc_exp <- read.table("data/exp.ACRG/ACRG-GC.txt")
meta <- read.csv("data/metaData/ACRG-Fig1.CSV")

# 2. Metadata pre-processing -----
# Mol. Subtype: 0=MSS/TP53-, 1=MSS/TP53+, 2 = MSI, 3= EMT
meta[which(meta$Mol.subtype == 0),"Mol.subtype"] <- "MSS/TP53-"
meta[which(meta$Mol.subtype == 1),"Mol.subtype"] <- "MSS/TP53+"
meta[which(meta$Mol.subtype == 2),"Mol.subtype"] <- "MSI"
meta[which(meta$Mol.subtype == 3),"Mol.subtype"] <- "EMT"

head(meta)

# 3. Merge datasets -----
library(dplyr)
gc_exp_t <- t(gc_exp) %>% as.data.frame()
gc_exp_t$Sample.ID <- rownames(gc_exp_t)
tmp <- unlist(str_split(gc_exp_t$Sample.ID, "_"))
gc_exp_t$Sample.ID <- c(tmp[seq(2,length(tmp),2)])

merged_gc <- merge(meta,gc_exp_t,by = "Sample.ID")
write.table(merged_gc, file = "data/exp.ACRG/ACRG-GC-merged-meta.txt")
