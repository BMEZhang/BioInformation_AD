# step3 结果可视化

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(stringr)

rm(list = ls()) #清空变量
gc()


### 1.数据读取------------------------------------------------------------------
output_dir = './DATA/'
gse_number = 'GSE83500'
filename <-  paste(output_dir, gse_number, '/step2output.Rdata', sep = "", collapse = NULL)
load(file = filename)

probe2entrez <- annot_table %>%
  dplyr::select(PROBEID = 'ID', ENTREZID = 'Entrez Gene') %>%
  filter(!is.na(ENTREZID) & ENTREZID != "")

# probe2entrez <- probe2entrez[!is.na(probe2entrez$PROBEID) & !is.na(probe2entrez$ENTREZID), ] 
multi_mapping_probes <- probe2entrez$PROBEID[grepl("///", probe2entrez$ENTREZID)]
probe_annot_filtered <- probe2entrez[!probe2entrez$PROBEID %in% multi_mapping_probes, ]
# 筛选 DEGs
deg_significant <- subset(all_deg_results, P.Value < 0.05 & abs(logFC) > 0)
# 查看筛选出的DEG数量
nrow(deg_significant)
# 查看上调和下调的基因数量
table(deg_significant$logFC > 0)
gene_list <- rownames(deg_significant) # 提取显著DEGs的基因名
converted <- probe_annot_filtered %>%
  filter(PROBEID %in% gene_list)
gene_entrez_ids <- unique(na.omit(converted$ENTREZID))
### 1.数据读取------------------------------------------------------------------



### 2.GO富集分析---------------------------------------------------------------
# GO富集分析（以 Biological Process 为例）
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
# 同样，可以分析MF和CC
# ego_mf <- enrichGO(... , ont = "MF")
# ego_cc <- enrichGO(... , ont = "CC")
### 2.GO富集分析---------------------------------------------------------------



### 3.KEGG通路------------------------------------------------------------------
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
### 3.KEGG通路------------------------------------------------------------------



### 4.可视化富集结果------------------------------------------------------------
# 4.1 条形图 (Barplot)
barplot(ego_bp, 
        showCategory = 15, # 显示最显著的前15个Term
        title = "GO Biological Process Enrichment",
        font.size = 8)

# 4.2 点图 (Dotplot) - 更常用，信息更丰富
dotplot(ego_bp, showCategory = 15)
dotplot(kk, showCategory = 15, title = "KEGG Pathway Enrichment")

# 4.3 网络图 (Enrichment Map) - 展示Term之间的重叠基因
ego_bp2 <- pairwise_termsim(ego_bp)
emapplot(ego_bp2, showCategory = 20)

# 4.4 基因-通路关系图 (Gene-Concept Network) - 需要Term数量较少
cnetplot(ego_bp, 
         showCategory = 5, # 只显示前5个Term，否则图会很乱
         circular = FALSE, 
         colorEdge = TRUE,
         node_label = "all") # 显示所有节点标签

# 4.5 山脊图 (Ridgeline Plot)
ridgeplot(ego_bp, showCategory = 15)
### 4.可视化富集结果------------------------------------------------------------



### 5.PPI分析-------------------------------------------------------------------
library(STRINGdb)
# 创建一个STRINGdb对象（针对人类物种）
string_db <- STRINGdb$new(version = "11.5", 
                          species = 9606, # 人的物种编号
                          score_threshold = 400, # 互作置信度阈值(0-1000)
                          input_directory = "")

# 映射我们的基因（使用之前转换好的ENTREZID数据框）
mapped_genes <- string_db$map(converted, "ENTREZID", removeUnmappedRows = TRUE)

# 获取互作对（边列表）
ppi_edges <- string_db$get_interactions(mapped_genes$STRING_id)

# 为边列表添加原始的Gene Symbol信息
# 首先创建一个映射字典
string_to_symbol <- setNames(mapped_genes$SYMBOL, mapped_genes$STRING_id)
ppi_edges$from_symbol <- string_to_symbol[ppi_edges$from]
ppi_edges$to_symbol <- string_to_symbol[ppi_edges$to]

# 保存PPI边列表，用于导入Cytoscape或其他软件
write.csv(ppi_edges, "PPI_Edge_List_for_Cytoscape.csv", row.names = FALSE)
### 5.PPI分析-------------------------------------------------------------------




# 保存结果------------------------------------------------------
write.csv(as.data.frame(ego_bp), "GO_BP_Enrichment_Results.csv", row.names = FALSE)
write.csv(as.data.frame(kk), "KEGG_Enrichment_Results.csv", row.names = FALSE)