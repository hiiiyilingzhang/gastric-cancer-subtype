# 炎症小体相关基因: NLRP3  AIM2  NLRC4  NLRP1  NLRP6  
# N3炎症小体：Caspase-1  GSDMD   ASC  IL-1β   IL-17A 
# 趋化、清道夫受体：CD36  CXCL20
# 中性粒细胞诱捕网 ：MPO  Ly6G5C  Ly6G5B  PAD4  PSGL1

# 1. Subset data -----
library(dplyr)
merged_gc <- read.table("data/exp.ACRG/ACRG-GC-merged-meta.txt",header = T)
subset.m <- select(merged_gc, one_of(c("Mol.subtype","NLRP3","AIM2","NLRC4","NLRP1","NLRP6","CASP1","GSDMD",
                                       "ASC","IL1B","IL17A","CD36","CXCL20","MPO","LY6G5C",
                                       "LY6G5B","PAD4","PSGL1")))
# Warning message:
#   Unknown columns: `ASC`, `CXCL20`, `LY6G5B`, `PAD4`, `PSGL1`

# 2. Create Metadata for heapmap -----
rownames(subset.m) <- paste("Patient",1:300,sep = "-")

annotation_col <- data.frame(ACRG.Subtype = subset.m$Mol.subtype)
rownames(annotation_col) <- paste("Patient",1:300,sep = "-")
head(annotation_col)

annotation_row = data.frame(
  Function = factor(rep(c("inflammasome-related genes", "N3 inflammasome", 
                          "Chemokine and scavenger receptor", "NETs"), 
                        c(5, 4, 1, 2)))
)
rownames(annotation_row) <- c("NLRP3","AIM2","NLRC4","NLRP1","NLRP6","CASP1","GSDMD","IL1B",
                              "IL17A","CD36","MPO","LY6G5C")
head(annotation_row)

# 3. Heatmap -----
library(pheatmap)
subset.m <- arrange(subset.m, Mol.subtype)
pdf(file= "results/ACRG-res.pdf")
pheatmap(t(subset.m[,-1]), color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_cols = FALSE, cluster_rows = FALSE, show_colnames=F,
         annotation_col = annotation_col, annotation_row = annotation_row)
dev.off()
