# 1. Read annoted data and metadata
exp_anno <- readRDS("data/exp.TCGA/TCGA-annoted.RDS")
metadata <- read.csv("data/metaData/TCGA-Fig1.CSV",header = T)

# 2. Rearrange Sample.ID and metadata
exp_anno_t <- t(exp_anno) %>% as.data.frame()
exp_anno_t$Sample.ID <- rownames(exp_anno_t)
# exp_anno_t[1:4,1:6]

exp_anno_t$Sample.ID <- str_sub(exp_anno_t$Sample.ID,start = 1,end=12)
exp_anno_t$Sample.ID <- gsub("\\.","\\-",exp_anno_t$Sample.ID)

metadata <- metadata[,-4]

# 3. Merge
merged_gc <- merge(metadata, exp_anno_t, by = "Sample.ID")
write.table(merged_gc,file = "data/exp.TCGA/TCGA-merged-metadata.txt",quote = F,row.names = F,sep = "\t")
