## GEO_DEGs_Analysis------------------------------------------------------------
# Integrating TSPO PET imaging and transcriptomics to unveil the role of neuroinflammation and amyloid‐β deposition in Alzheimer’s disease
# PMID: 37801139

rm(list = ls()) #清空变量
gc()
source('function_zwr.R')



### 0.数据输入------------------------------------------------------------------
library_zwr()
gse_number = 'GSE83500'
output_dir1 <- paste('D:/Bio_Data/GEO/', gse_number, '/',sep = "", collapse = NULL)

GROUP_HC <- "HC"
GROUP_EXP <- "EXP"
GROUP_criterion <- "non-MI"
# ------------------------------------------------------------------------------


### 1.数据读取------------------------------------------------------------------
data_temp <- data_download_zwr(output_dir1,gse_number) # 下载数据，并做标准化检查
sample_matrix <- data_temp[[1]]
sample_metadata <- data_temp[[2]]

gpl_number <- data_temp[[3]]@annotation # 如果是GPL列表，取第一个
gpl <- getGEO(gpl_number, destdir = output_dir1)  # 下载平台信息
annot_table <- Table(gpl)   # 提取原始注释表

probe2entrez <- annot_table %>%
  dplyr::select(PROBEID = 'ID', SYMBOL = 'Gene Symbol', ENTREZID = 'Entrez Gene') %>%
  filter(!is.na(ENTREZID) & ENTREZID != "")
multi_mapping_probes <- probe2entrez$SYMBOL[grepl("///", probe2entrez$ENTREZID)]
probe_annot_filtered <- probe2entrez[!probe2entrez$SYMBOL %in% multi_mapping_probes, ]
probe_annot_filtered <- probe_annot_filtered %>%
  mutate(
    ENTREZID = na_if(ENTREZID, ""),     # 将空字符串 "" 转换为 NA
    ENTREZID = na_if(ENTREZID, "---"),  # 将 "---" 转换为 NA
    ENTREZID = na_if(ENTREZID, "unknown") # 如果需要，还可以处理其他标记
  ) %>%
  na.omit() # 现在，na.omit() 会一次性删除所有包含NA的行
examine_gpl_zwr(probe_annot_filtered) # 检查GPL去空后数据
expr_mat <- desert_expr_zwr(sample_matrix, probe2entrez) # 舍弃多映射探针的表达矩阵

metadata_Name <- paste(output_dir1, gse_number, "_metadata.csv", sep = "", collapse = NULL)
write.csv(sample_metadata, metadata_Name, row.names = F)
exprdata_Name <- paste(output_dir1, gse_number, "_matrix_expr.csv", sep = "", collapse = NULL)
write.csv(expr_mat, exprdata_Name, row.names = F)
filename <-  paste(output_dir1, 'step1output.Rdata', sep = "", collapse = NULL)
save(gse_number,sample_metadata,annot_table,expr_mat,gpl_number,file = filename) # 保存数据
# ------------------------------------------------------------------------------


### 2.数据分组处理----------------------------------------------------------------
SYMBOL_ids <- rownames(expr_mat) 
Group=ifelse(str_detect(sample_metadata$title,GROUP_criterion),
             GROUP_HC,
             GROUP_EXP)
Group_ID <- paste0(Group, ave(Group, Group, FUN = seq_along))
Group_mat <- data.frame(sample_metadata$geo_accession,Group_ID = Group_ID,Group = Group)
colnames(Group_mat) <- c('GSM', 'SampleID', 'Group')
colnames(expr_mat) <- Group_mat$SampleID
examine_meta_zwr(expr_mat, Group_mat) # 检查表达矩阵和样本信息是否匹配
# ------------------------------------------------------------------------------


### 3.拟合线性模型--------------------------------------------------------------
all_deg_results <- fit_line_zwr(Group_mat,expr_mat)
# 导出所有基因的完整结果
file_All_Gene <- paste(output_dir1, gse_number, 'All_Gene_Results_EXP_vs_HC.csv', sep = "", collapse = NULL)
write.csv(all_deg_results, file = file_All_Gene, row.names = TRUE)
# ------------------------------------------------------------------------------


### 4.搜索差异基因--------------------------------------------------------------

# 筛选 DEGs
deg_significant <- subset(all_deg_results, P.Value < 0.05 & abs(logFC) > 0)

# # 查看筛选出的DEG数量
# nrow(deg_significant)
# 查看上调和下调的基因数量
table(deg_significant$logFC > 0)
deg_gene_list <- rownames(deg_significant) # 提取显著DEGs的基因名
converted <- probe_annot_filtered %>%
  filter(SYMBOL %in% deg_gene_list)
SYMBOL_ids <- unique(na.omit(converted$SYMBOL))
gene_entrez_ids <- unique(na.omit(converted$ENTREZID))
deg_significant_gene <- data.frame(SYMBOL = SYMBOL_ids, ENTREZID = gene_entrez_ids)
# 导出显著的DEGs GENE列表
deg_sig_Name <- paste(output_dir1, gse_number, "_DEG_SIG.csv", sep = "", collapse = NULL)
write.csv(deg_significant_gene, deg_sig_Name, row.names = F)
# 导出显著的DEGs列表
file_DEGs <- paste(output_dir1, gse_number, 'Significant_DEGs_EXP_vs_HC_Pval_0.05_logFC_0.csv', sep = "", collapse = NULL)
write.csv(deg_significant, file = file_DEGs, row.names = TRUE)
# ------------------------------------------------------------------------------


### 获取交集基因进行后续分析(Venn plot)-----------------------------------------
deg_list1 <- rownames(deg_significant[30:180,]) # 例如 EXP vs HC 的全部DEGs
deg_list2 <- rownames(deg_significant[125:200,]) # 例如 Male vs Female 的DEGs
deg_list3 <- rownames(deg_significant[161:250,]) # 例如 Early AD vs Late AD 的DEGs
overlaps <- venn(list(ADvsHC = deg_list1, MvsF = deg_list2, EarlyvsLate = deg_list3), show.plot = FALSE)
# 获取特定交集的基因名
intersections <- attr(overlaps, "intersections")
# 查看所有交集的名称
names(intersections)
# 获取特定交集的基因
common_Symbol <- intersections$`ADvsHC:MvsF`
# 获取只在某一个列表中出现的基因
only_in_ad_hc <- setdiff(deg_list1, union(deg_list2, deg_list3))
# ------------------------------------------------------------------------------


### 5.GO富集分析----------------------------------------------------------------
# GO富集分析
ego_bp <- enrichGO(gene = gene_entrez_ids,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", # 输入的基因ID类型
                   ont = "BP", # "BP", "MF", "CC"
                   pAdjustMethod = "BH", # 多重检验校正方法
                   pvalueCutoff = 0.05, # P值阈值
                   qvalueCutoff = 0.2, # Q值阈值
                   readable = TRUE) # 是否将ENTREZID转换为Symbol，方便解读结果
# 查看结果摘要
head(ego_bp, 10)
filename <- paste(output_dir1, 'GO_BP_Enrichment_Results.csv', sep = "", collapse = NULL)
write.csv(as.data.frame(ego_bp), filename, row.names = FALSE)
# ------------------------------------------------------------------------------


### 6.KEGG通路------------------------------------------------------------------
# KEGG富集分析
kk <- enrichKEGG(gene = gene_entrez_ids,
                 organism = "hsa", # 物种缩写：hsa（人），mmu（小鼠）等
                 keyType = "kegg", 
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)

# 将ENTREZID转换为Gene Symbol，便于解读
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
# 查看结果摘要
head(kk, 10)
filename <- paste(output_dir1, 'KEGG_Enrichment_Results.csv', sep = "", collapse = NULL)
write.csv(as.data.frame(kk), filename, row.names = FALSE)
# ------------------------------------------------------------------------------


### 7.可视化富集结果------------------------------------------------------------
#可视化选项 1 = "barplot"; 2 = "dotplot"; 3 = "ego_bp2"; 4 = "cnetplot"
vis_opt = 2
vis_GOKEGG_zwr(vis_opt,ego_bp,kk)
# ------------------------------------------------------------------------------


### 8.PPI分析-------------------------------------------------------------------
ppi_edges <- PPI_analysis_zwr(deg_significant_gene)
# ensembl_ids <- mapIds(org.Hs.eg.db,
#                       keys = gene_entrez_ids,
#                       column = "ENSEMBL",
#                       keytype = "ENTREZID",
#                       multiVals = "first")  # 处理一个ENTREZID对应多个ENSEMBL的情况
# ppi_edges <- PPI_analysis_zwr(gene_entrez_ids)
# ppi_edges2 <- PPI_analysis_zwr(ensembl_ids)
ppi_name <- paste(output_dir1, 'PPI_Edge_List_for_Cytoscape.csv', sep = "", collapse = NULL)
write.csv(ppi_edges, ppi_name, row.names = FALSE)
# ------------------------------------------------------------------------------


### 9.绘制火山图 (Volcano Plot)-------------------------------------------------
Volcano_zwr(all_deg_results,P_value=0.05,logFC_value=0)
# ------------------------------------------------------------------------------


### 10.绘制热图 (Heatmap)--------------------------------------------------------
HeatMap_zwr(deg_significant,expr_mat, Group_mat)
# ------------------------------------------------------------------------------


### 11.韦恩图 (Venn Diagram)-----------------------------------------------------
list_deg <- list(deg_significant)
Venn_zwr(list_deg)
# ------------------------------------------------------------------------------


### 12.绘制UpSet图（对于比较多个集合）------------------------------------------
gene_lists <- list(
  ADvsHC = deg_list1,
  MvsF = deg_list2,
  EarlyvsLate = deg_list3
)
# 创建输入数据
input_data <- fromList(gene_lists)
# 绘制UpSet图
upset(input_data, 
      sets = c("ADvsHC", "MvsF", "EarlyvsLate"),
      order.by = "freq", 
      mainbar.y.label = "Number of Common Genes",
      sets.x.label = "Genes per List")
# ------------------------------------------------------------------------------


### 13.绘制MA图 (MA Plot)-------------------------------------------------------
Ma_zwr(Group_mat, expr_mat)
# ------------------------------------------------------------------------------
