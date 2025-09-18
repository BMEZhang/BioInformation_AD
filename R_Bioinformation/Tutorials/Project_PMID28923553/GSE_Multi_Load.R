# 多数据处理
## GEO_DEGs_Analysis------------------------------------------------------------
# A systematic integrated analysis of brain expression profiles reveals YAP1 and other prioritized hub genes as important upstream regulators in Alzheimer's disease
# PMID: 28923553
# 2025/09/18

rm(list = ls()) #清空变量
gc()
source('function_zwr.R')

### 0.Input---------------------------------------------------------------------
gse_number = 'GSE66333'
output_dir1 <- paste('D:/Bio_Data/GEO/DATA_PMID28923553/', gse_number, '/',sep = "", collapse = NULL)
filename <-  paste(output_dir1, gse_number, '_step1output.Rdata', sep = "", collapse = NULL)
load(filename)
# ------------------------------------------------------------------------------





