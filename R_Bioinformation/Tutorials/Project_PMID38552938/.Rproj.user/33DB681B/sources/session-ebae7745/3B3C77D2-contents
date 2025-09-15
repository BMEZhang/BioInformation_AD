# step 0
# get data
# ("myocardial infarction"[MeSH Terms]) AND ("Homo sapiens"[porgn])
# AND ("gse"[Filter]) AND ("Expression profiling by array"[Filter])
# search 159 items

rm(list = ls()) #清空变量
gc()
library(GEOquery)
library(Biobase)
# data <- read.table("GSE83500_series_matrix.txt.gz", header = TRUE, sep = "\t", comment.char = "!", row.names = 1)
# 37个样本，49386个基因

###########################################################
#检查数据是否已经标准化
examine_matrix <- function(exp) {
  #exp = log2(exp+1) #看数据是否已经取过log值——数值是否跨数量级
  ex <- exp #把需要判断的矩阵命名为ex
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) { ex[which(ex <= 0)] <- NaN
  exp <- log2(ex+1)
  print("log2 transform finished")}else{print("log2 transform not needed")}
  boxplot(exp) #看一下数据是否基本在一条直线上
  return(exp)
}

###########################################################

### Download GEO database
output_dir = './DATA/'
gse_number = 'GSE83500'
output_dir <- paste(output_dir, gse_number, '/',sep = "", collapse = NULL)
dir.create(output_dir)
eSet <- getGEO(gse_number, 
               destdir = output_dir, 
               AnnotGPL = F,
               getGPL = F)
length(eSet) #看有几个GSE
SetName <- paste(output_dir, gse_number, '.gset.Rdata', sep = "", collapse = NULL)
save(eSet,file = SetName)

# 1.提取表达矩阵
sample_matrix <- exprs(eSet[[1]])
head(sample_matrix)[, 1:5]
sample_matrix <- examine_matrix(sample_matrix) #检查数据是否已经标准化


# 2.提取样本元数据
sample_metadata <- pData(eSet[[1]])

# 3.调整metadata的行名顺序与exp列名完全一致
p = identical(rownames(sample_metadata),colnames(sample_matrix));p
if(!p) sample_matrix = sample_matrix[,match(rownames(sample_metadata),colnames(sample_matrix))]
metadata_Name <- paste(output_dir, "GSE83500_metadata.csv", sep = "", collapse = NULL)
write.csv(sample_metadata, metadata_Name, row.names = F)

# 4.提取芯片平台编号
gpl_number <- eSet[[1]]@annotation # 如果是GPL列表，取第一个
filename <-  paste(output_dir, 'step0output.Rdata', sep = "", collapse = NULL)
save(gse_number,sample_metadata,sample_matrix,gpl_number,file = filename)


