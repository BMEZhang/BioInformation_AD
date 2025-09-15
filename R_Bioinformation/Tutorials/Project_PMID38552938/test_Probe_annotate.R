# # test_Probe_annotate
# 
# output_dir = './DATA/'
# gse_number = 'GSE83500'
# filename <-  paste(output_dir, gse_number, '/step1output.Rdata', sep = "", collapse = NULL)
# load(file = filename)
# 
# 
# options("repos"="https://mirrors.ustc.edu.cn/CRAN/")
# if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)
# options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("hgu133plus2.db")
# library(hgu133plus2.db)
# 
# probe_ids <- rownames(sample_matrix)  # 假设exprs_data是你的表达矩阵，行名为探针ID
# annotations <- select(hgu133plus2.db, keys=probe_ids, columns=c("SYMBOL", "ENTREZID", "GENENAME"), keytype="PROBEID")
# 
# 
# 
# 
# library(GEOquery)
# gpl <- getGEO("GPL570", destdir = ".")  # 下载平台信息
# annot_table <- Table(gpl)   # 提取注释表
# colnames(annot_table)       # 查看可用列
# probe_annot <- annot_table[, c("ID", "Gene Symbol", "ENTREZ_GENE_ID")]  # 选择需要的列（通常是ID和Gene Symbol）列名可能因平台而异
# colnames(probe_annot) <- c("PROBEID", "SYMBOL", "ENTREZID")  # 标准化列名
# probe_annot <- probe_annot[probe_annot$SYMBOL != "" & !is.na(probe_annot$SYMBOL), ]   # 去除空值或无效注释




# head(keys(hgu133plus2.db, keytype = "SYMBOL"), 10)  # 前5个探针
# [1] "A1BG"     "A2M"      "A2MP1"    "NAT1"     "NAT2"       
# > head(keys(hgu133plus2.db, keytype = "ENTREZID"), 10)  # 前5个探针
# [1] "1"         "10"        "100"       "1000"      "10000" 
# > head(keys(hgu133plus2.db, keytype = "GENENAME"), 10)  # 前5个探针
# [1] "alpha-1-B glycoprotein"                  "alpha-2-macroglobulin"                   "alpha-2-macroglobulin pseudogene 1"     
# [4] "N-acetyltransferase 1"                   "N-acetyltransferase 2"
# head(keys(hgu133plus2.db, keytype = "PROBEID"), 10)  # 前5个探针
# [1] "1007_s_at" "1053_at"   "117_at"    "121_at"    "1255_g_at"


rm(list = ls())


# 加载示例数据（假设expr_matrix为表达矩阵，行名为探针ID，列名为样本）
expr_matrix <- matrix(rnorm(1000), nrow=100, dimnames=list(paste0("Probe", 1:100), paste0("Sample", 1:10)))

# 创建探针与基因的映射表（示例：随机生成100个探针对应50个基因）
set.seed(123)
gene_symbols <- paste0("Gene", 1:50)
probe2gene <- data.frame(
  probe = rownames(expr_matrix),
  gene = sample(gene_symbols, 100, replace = TRUE)
)


# 安装并加载WGCNA
if (!require("WGCNA")) install.packages("WGCNA")
library(WGCNA)

# 假设expr_matrix为表达矩阵，probe2gene为探针-基因映射表
collapse_result <- collapseRows(
  datET = expr_matrix,
  rowGroup = probe2gene$gene,      # 基因名向量（与表达矩阵行顺序一致）
  rowID = rownames(expr_matrix)    # 探针名向量
)

# 提取合并后的矩阵
expr_collapsed <- collapse_result$datETcollapsed