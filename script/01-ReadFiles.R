# 1. Read files -----
library(affy)
dir_cels <- 'data/rawData/untarData/'   # path to your .cel data
affy_data <- ReadAffy(celfile.path=dir_cels)

# 2. Set rownames -----
library(stringr)

## extract sample name
raw.names<-sampleNames(affy_data)

## separated based on '_' and extract 10th characters
## e.g. 'GSM1523838_INT201201-10-1_443T_HG-U133_Plus_2_.CEL.gz' can be separated by '_' as:
## GSM1523838, INT201201-10-1, 443T,HG-U133, Plus, 2, .CEL.gz
new<-str_split_fixed(raw.names, "_", 9)

new.1<- new[1:150,]
new.1 <- new.1[,-9]    # remove V9
new.1[,3] <- str_remove(string = new.1[,3],pattern = "T")    # remove "T" ended in V3
new.1[78,3] <- 441    # '441-2' --> '441'
# head(new.1)

new.2 <- new[151:300,]
head(new.2)    # 8 columns and V3 'tr' need to be removed
new.2[1,3:8] <- c("Tr","207","HG-U133","Plus","2",".CEL.gz")
new.2 <- new.2[,-3]    # remove V3, need to consistent with new.1
# head(new.2)

new <- as.data.frame(rbind(new.1,new.2))
new <- tidyr::unite(new, col = "sampleName", V1, V3,sep = "_")
head(new)
new.names<-new[,1]

sampleNames(affy_data)<-new.names
pData(affy_data)    # check data

# 3. Data Normalization -----
est.mas5 <- mas5(affy_data)

# 4. Extract expression data -----
gc_exp <- exprs(est.mas5)
gc_exp[1:4,1:5]

# 5. Probe annotation -----
# BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)

ids <- toTable(hgu133plus2SYMBOL)
head(ids)

table(rownames(gc_exp) %in% ids$probe_id)
# FALSE  TRUE 
# 11526 43149 
# 43149 probe can be annotated

data <- gc_exp[rownames(gc_exp) %in% ids$probe_id,]    # 取出原始数据中有注释的探针
dim(data)

ids <- ids[match(rownames(data),ids$probe_id),]    # 取出要整合的探针id
tmp <- by(data,ids$symbol,function(x) rownames(x)[which.max(rowMeans(x))] )    # 去重，如果一个探针号对应多个基因，就取表达值最高的
probes<-as.character(tmp)
data<-data[rownames(data) %in% probes,]
ids<-ids[ids$probe_id %in% probes,]
row.names(data)<-ids$symbol

head(data)
dim(data)

# 6. Save data -----
write.table(data, file="data/exp.ACRG/ACRG-GC.txt") 
