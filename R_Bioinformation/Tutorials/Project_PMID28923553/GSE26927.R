## GEO_DEGs_Analysis------------------------------------------------------------
# A systematic integrated analysis of brain expression profiles reveals YAP1 and other prioritized hub genes as important upstream regulators in Alzheimer's disease
# PMID: 28923553
# 2025/09/17

rm(list = ls()) #清空变量
gc()
source('function_zwr.R')


### 0.数据输入------------------------------------------------------------------
gse_number = 'GSE33000'
output_dir1 <- paste('D:/Bio_Data/GEO/DATA_PMID28923553/', gse_number, '/',sep = "", collapse = NULL)
# ------------------------------------------------------------------------------


### 1.数据读取------------------------------------------------------------------
data_temp <- data_download_zwr(output_dir1,gse_number) # 下载数据，并做标准化检查
dev.off()
# 提取表达矩阵
sample_matrix <- data_temp[[1]] 
# 提取元数据矩阵
sample_metadata <- data_temp[[2]] # pdata
# 提取原始注释表
annot_table <- data_temp[[3]]  
rm(data_temp)
# 提取PROBEID、SYMBOL以及ENTREZID注释表
colnames(annot_table)
probe_annot <- probe_annot_zwr(annot_table,"ID","Gene Symbol","ENTREZ_GENE_ID" )
cat(sprintf("注释表包含: %d 个独立SYMBOL\n", length(unique(probe_annot$SYMBOL))))
# 按探针注释表提取
expr_mat <- desert_expr_zwr(sample_matrix, probe_annot) # 舍弃多映射探针的表达矩阵以及注释矩阵
cat(sprintf("表达矩阵包含: %d 个SYMBOL, 和 %d 个Sample\n",dim(expr_mat)[1], dim(expr_mat)[2]))
# 保存expr_mat，sample_metadata，probe_annot至csv表以及.Rdata文件。
save_mat_zwr(output_dir1, gse_number,expr_mat, sample_metadata)
# ------------------------------------------------------------------------------


### 2.数据分组处理----------------------------------------------------------------
# 第1步：获取数据并创建因子
pdata <- sample_metadata
# 自定义分组，与数据相关
pdata <- pdata %>%
  mutate(
    Disease = case_when(
      grepl("AD", title) ~ "AD",
      grepl("ALS", title) ~ "ALS",
      grepl("HD", title) ~ "HD",
      grepl("MS", title) ~ "MS",
      grepl("PD", title) ~ "PD",
      grepl("SCHIZ", title) ~ "SCHIZ",
      TRUE ~ NA_character_
    ),
    Status = case_when(
      grepl("Control", title) ~ "Control",
      grepl("Disease", title) ~ "Disease",
      TRUE ~ NA_character_
    )
  ) %>%
  # 确保没有NA值
  filter(!is.na(Disease), !is.na(Status)) 

pdata$Group <- interaction(pdata$Disease, pdata$Status, sep = "_")
# 将因子转换为因子类型并设置参考水平（通常对照组或时间0点为参考）
pdata$Disease <- factor(pdata$Disease, levels = c("AD", "ALS", "HD", "MS", "PD", "SCHIZ"))
pdata$Status <- factor(pdata$Status, levels = c("Control", "Disease"))
pdata$Group <- factor(pdata$Group)
# 查看分组表，确认样本数量
table(pdata$Disease, pdata$Status)

# 第2步：构建设计矩阵 (Design Matrix)
# 方法A：使用交互分组因子（最直接、最灵活）
design <- model.matrix(~ 0 + pdata$Group)
colnames(design) <- levels(pdata$Group) # 重命名列名为分组名
rownames(design) <- rownames(pdata) # 重命名行名为样本名
print(design)

expr_mat <- na.omit(expr_mat)
examine_meta_zwr(expr_mat, pdata) # 检查表达矩阵和样本信息是否匹配
dev.off()

# 第3步：定义对比矩阵 (Contrast Matrix)
# 基于方法A的设计矩阵
fit <- lmFit(expr_mat, design)
# 定义对比
contrast.matrix <- makeContrasts(
  # 1. 简单效应：AD
  AD_Disease_vs_Control = AD_Disease - AD_Control,
  # 2. 简单效应：ALS
  ALS_Disease_vs_Control = ALS_Disease - ALS_Control,
  # 3. 简单效应：HD
  HD_Disease_vs_Control = HD_Disease - HD_Control,
  # 4. 简单效应：MS
  MS_Disease_vs_Control = MS_Disease - MS_Control,
  # 5. 简单效应：PD
  PD_Disease_vs_Control = PD_Disease - PD_Control,
  # 6. 简单效应：SCHIZ
  SCHIZ_Disease_vs_Control = SCHIZ_Disease - SCHIZ_Control,
  
  # 7. 主效应: Disease vs Control (可以看作是所有时间点的平均效应)
  # (KO_0h + KO_12h + KO_24h)/3 - (WT_0h + WT_12h + WT_24h)/3
  Disease_vs_Control = (AD_Disease + ALS_Disease + HD_Disease + MS_Disease + PD_Disease + SCHIZ_Disease)/6 - (AD_Control + ALS_Control + HD_Control + MS_Control + PD_Control + SCHIZ_Control)/6,
  # 8. 交互效应: (KO在12h的效果 - KO在0h的效果) vs (WT在12h的效果 - WT在0h的效果)
  # 即: (KO_12h - KO_0h) - (WT_12h - WT_0h)
  Interaction_AD_ALS = (AD_Disease - AD_Control) - (ALS_Disease - ALS_Control),
  # 9. 交互效应: 同上，针对24h
  Interaction_HD_MS = (HD_Disease - HD_Control) - (MS_Disease - MS_Control),
  levels = design
)

# 第4步：拟合模型并计算结果
# 将对比矩阵代入拟合模型
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
# # 提取所有基因的统计结果
# all_deg_results <- topTable(fit2, number = Inf, coef = "EXPvsHC", adjust.method = "BH", sort.by = "P")


# 查看任意对比的结果
de_genes_AD <- topTable(fit2, coef = "AD_Disease_vs_Control", number = Inf, adjust.method = "BH", sort.by = "P")
de_genes_Interaction <- topTable(fit2, coef = "Disease_vs_Control", number = Inf, adjust.method = "BH", sort.by = "P")

# 对所有比较进行F检验（查看哪些基因在任何对比中差异表达）
results <- decideTests(fit2)
summary(results)

# ------------------------------------------------------------------------------


### 3.拟合线性模型--------------------------------------------------------------
# all_deg_results <- fit_line_zwr(Group_mat,expr_mat)
# 导出AD_vs_Control基因的完整结果
file_AD_Gene <- paste(output_dir1, gse_number, '_All_Gene_Results_AD_vs_HC.csv', sep = "", collapse = NULL)
write.csv(de_genes_AD, file = file_AD_Gene, row.names = TRUE)
# 导出Disease_vs_Control基因的完整结果
file_Disease_Gene <- paste(output_dir1, gse_number, '_All_Gene_Results_Disease_vs_HC.csv', sep = "", collapse = NULL)
write.csv(de_genes_Interaction, file = file_Disease_Gene, row.names = TRUE)

# ------------------------------------------------------------------------------


### 4.搜索差异基因--------------------------------------------------------------

# 筛选 DEGs
deg_significant_AD <- subset(de_genes_AD, P.Value < 0.05 & abs(logFC) > 0)
deg_significant_Disease <- subset(de_genes_Interaction, P.Value < 0.05 & abs(logFC) > 0)

# # 查看筛选出的DEG数量
# nrow(deg_significant)
# 查看上调和下调的基因数量
table(deg_significant_AD$logFC > 0)
deg_gene_AD_list <- rownames(deg_significant_AD) # 提取显著DEGs的基因名

converted <- probe_annot %>%
  filter(SYMBOL %in% deg_gene_AD_list)
SYMBOL_ids <- unique(na.omit(converted$SYMBOL))
gene_entrez_ids <- unique(na.omit(converted$ENTREZID))
deg_significant_AD_gene <- data.frame(SYMBOL = SYMBOL_ids, ENTREZID = gene_entrez_ids)
# 导出显著的DEGs GENE列表
deg_sig_Name <- paste(output_dir1, gse_number, "_DEG_SIG.csv", sep = "", collapse = NULL)
write.csv(deg_significant_AD_gene, deg_sig_Name, row.names = F)
# 导出显著的DEGs列表
file_DEGs <- paste(output_dir1, gse_number, '_Significant_DEGs_EXP_vs_HC_Pval_0.05_logFC_0.csv', sep = "", collapse = NULL)
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
filename <- paste(output_dir1, gse_number, '_GO_BP_Enrichment_Results.csv', sep = "", collapse = NULL)
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
filename <- paste(output_dir1, gse_number, '_KEGG_Enrichment_Results.csv', sep = "", collapse = NULL)
write.csv(as.data.frame(kk), filename, row.names = FALSE)
# ------------------------------------------------------------------------------


### 7.可视化富集结果------------------------------------------------------------
#可视化选项 1 = "barplot"; 2 = "dotplot"; 3 = "ego_bp2"; 4 = "cnetplot"
vis_opt = 4
vis_GOKEGG_zwr(vis_opt,ego_bp,kk,output_dir1,gse_number)
# ------------------------------------------------------------------------------


### 8.PPI分析-------------------------------------------------------------------
ppi_edges <- PPI_analysis_zwr(deg_significant_AD_gene)
# ensembl_ids <- mapIds(org.Hs.eg.db,
#                       keys = gene_entrez_ids,
#                       column = "ENSEMBL",
#                       keytype = "ENTREZID",
#                       multiVals = "first")  # 处理一个ENTREZID对应多个ENSEMBL的情况
# ppi_edges <- PPI_analysis_zwr(gene_entrez_ids)
# ppi_edges2 <- PPI_analysis_zwr(ensembl_ids)
ppi_name <- paste(output_dir1, gse_number, '_PPI_Edge_List_for_Cytoscape.csv', sep = "", collapse = NULL)
write.csv(ppi_edges, ppi_name, row.names = FALSE)
# ------------------------------------------------------------------------------


### 9.绘制火山图 (Volcano Plot)-------------------------------------------------
Volcano_zwr(de_genes_AD,P_value=0.05,logFC_value=0)
# ------------------------------------------------------------------------------


### 10.绘制热图 (Heatmap)--------------------------------------------------------
HeatMap_zwr(de_genes_AD,expr_mat, pdata)
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
