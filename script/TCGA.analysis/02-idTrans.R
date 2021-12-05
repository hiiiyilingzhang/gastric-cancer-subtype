library(dplyr)
library(tidyr)
library(rtracklayer)

# # TCGA的基因组由于是之前注释的，所以是gencode v22版本
# # GTF的下载从官方网站，https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz
# gtf1 <- rtracklayer::import('data/Reference/Homo_sapiens.GRCh38.100.gtf')
# gtf_df <- as.data.frame(gtf1)
# 
# # 选择两列，有2908224行
# gtf_df1<-select(gtf_df,gene_id,gene_name)
# # 去除重复行，有60683行，但里面还有多个ID对应一个基因的行
# gtf_df<-distinct(gtf_df1,gene_id,gene_name)
# 
# save(gtf_df,file="data/Reference/GTF_selected.Rdata")
# rm(gtf_df1,gtf_df2,gtf_df,gtf1)
# 
# # 进行注释替换
# mRNAmatrix<-read.table("data/2.ExprData/mRNAmatrix.txt",sep = "\t",stringsAsFactors = F,header = T)
# load("../../CancerGenomics/TCGA/tcga-classproj/data/Reference/GTF_selected.Rdata")
# 
# # 表达矩阵中的ensemble号中含有小数点， 而GTF文件中没有,调整一下，先以“.”分列，在去掉小数点后的列
# # 此行看下名字，不一定需要colnames(mRNAmatrix)<-”gene_id” 
# mRNAmatrix_nopoint <- mRNAmatrix %>% tidyr::separate(gene_id,into = c("gene_id","drop"),sep="\\.") %>% dplyr::select(-drop)
# 
# # 采用左连接，左边的全部保留
# mRNA_symbol<-left_join(x=mRNAmatrix_nopoint,y=gtf_df,by="gene_id")
# # 删除没能注释的行
# mRNA_symbol<-na.omit(mRNA_symbol)
# # 构建新的表达矩阵，第一列为gene symbol
# mRNA_symbol2<-cbind(gene_name=mRNA_symbol$gene_name,mRNA_symbol[,-c(1,ncol(mRNA_symbol))])
# 
# save(mRNA_symbol2,file="data/3.Symbol/mRNA_symbol.Rdata")
# write.table(mRNA_symbol2,file = "data/3.Symbol/symbol.txt",quote = F,row.names = F,sep = "\t")

# -----------------------------------------------------------------------

mRNAmatrix<-read.table("data/exp.TCGA/mRNAmatrix.txt",sep = "\t",stringsAsFactors = F,header = T)

# https://github.com/xjsun1221/tinyarray
library(tinyarray)  

rownames(mRNAmatrix) <- mRNAmatrix$gene_id
mRNAmatrix <- select(mRNAmatrix, -gene_id)
exp_anno <- trans_exp(mRNAmatrix,mrna_only = T)

saveRDS(exp_anno,file="data/exp.TCGA/TCGA-annoted.RDS")
write.table(exp_anno,file = "data/exp.TCGA/TCGA-annoted.txt",quote = F,row.names = F,sep = "\t")
