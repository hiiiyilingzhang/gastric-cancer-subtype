#解析meta文件，目的是为了拿到 样本名-barcode对应关系,
#注意json文件的位置
metadata <- jsonlite::fromJSON("data/rawData.TCGA/metadata.cart.2021-12-05.json")

library(dplyr)

# 转换的信息是两列filename和associate dentities
# 样本名称在associated_entities,列表中里面包括了 entity_id, case_id, entity_submitter_id, entity_type这四个项目
metadata_id <- metadata %>% dplyr::select(c(file_name,associated_entities)) 
# metadata$associatedentities[1][[1]]

# 把filename和 associated_entities中的entity_submitter_id提取出来，做成数据框批量对应转换
# substr(s, first, last)
# nchar(metadata_id$file_name[i])-3 去掉.gz
naid_df <- data.frame()
for (i in 1:dim(metadata_id)[1]){
  naid_df[i,1] <- substr(metadata_id$file_name[i],1,nchar(metadata_id$file_name[i])-3)
  naid_df[i,2] <- metadata_id$associated_entities[i][[1]]$entity_submitter_id
}

# 提取前12个
library(stringi)
library(stringr)
for (i in 1:dim(naid_df)[1]){
  naid_df$barcode[i]<- str_sub(naid_df[i,2],start=1,end=12)
}

# 提取出文件夹中的txt.gz
# 这里面有个技巧，先读文件夹中的内容，再新建文件夹“files”,要不然会报错
dir1 <- dir("data/rawData.TCGA/files/")
dir.create("data/rawData.TCGA/purefiles/") 
for (dirname in dir1){  
  file <- list.files(paste0("data/rawData.TCGA/files/",dirname),pattern = "*.FPKM.txt.gz")  #找到对应文件夹中的内容，pattern可以是正则表达式
  file.copy(paste0("data/rawData.TCGA/files/",dirname,"/",file),"data/rawData.TCGA/purefiles/")  #复制内容到新的文件夹
}
rm(dir1,dirname)

## gunzip
## remove MANIFEST.txt first
files1<-dir("data/rawData.TCGA/purefiles/")
library(R.utils)

for (i in 1:length(files1)){
  gunzip(paste0("data/rawData.TCGA/purefiles/",files1[i]),remove=TRUE)  ## 解压时原来的文件删掉
}
rm(files1,file,i)

#合并
library(data.table)

nameList <- list.files("data/rawData.TCGA/purefiles/",pattern = "*.FPKM.txt") ## 注意pattern 参数，获取该目录下所有的结尾为FPKM.txt
location <- which(naid_df==nameList[1],arr.ind = TRUE)  ## which函数有一个已知value返回坐标的功能
TCGA_id <- as.character(naid_df[location[1],2]) ## 通过坐标，获取TCGA_id

expr_df<- fread(paste0("data/rawData.TCGA/purefiles/",nameList[1])) ## 读入第一个文件，保存为data.frame
names(expr_df) <- c("gene_id",TCGA_id)  ## 给刚才数据库命名

for (i in 2:length(nameList)){
  location <- which(naid_df==nameList[i],arr.ind = TRUE)
  TCGA_id <- as.character(naid_df[location[1],2])
  dfnew <- fread(paste0("data/rawData.TCGA/purefiles/",nameList[i]))
  names(dfnew) <- c("gene_id",TCGA_id)
  expr_df <- inner_join(expr_df,dfnew,by="gene_id")
}

mRNAmatrix <-expr_df
saveRDS(mRNAmatrix,file = "data/exp.TCGA/mRNAmatrix.RDS")
write.table(mRNAmatrix,file = "data/exp.TCGA/mRNAmatrix.txt",quote = F,row.names = F,sep = "\t")
