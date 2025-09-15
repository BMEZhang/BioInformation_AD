# step3 结果可视化
library(limma)
library(ggplot2)
library(gplots)
library(pheatmap)
library(dplyr)
library(grid)
library(VennDiagram)
library(UpSetR)
library(ggrepel) # 防止标签重叠

rm(list = ls()) #清空变量
gc()



### 1.数据读取------------------------------------------------------------------
output_dir = './DATA/'
gse_number = 'GSE83500'
filename <-  paste(output_dir, gse_number, '/step2output.Rdata', sep = "", collapse = NULL)
load(file = filename)



### 2.应用阈值筛选差异表达基因--------------------------------------------------
# 筛选 DEGs
deg_significant <- subset(all_deg_results, P.Value < 0.05 & abs(logFC) > 0)
# 查看筛选出的DEG数量
nrow(deg_significant)
# 查看上调和下调的基因数量
table(deg_significant$logFC > 0)



### 3.结果可视化与导出----------------------------------------------------------

# 3.1 火山图 (Volcano Plot)------------------------------
# 为绘图数据添加一列，标记显著基因
all_deg_results$Significant <- ifelse(all_deg_results$P.Value < 0.05 & abs(all_deg_results$logFC) > 0, "Significant", "Not Sig")
# 为数据框添加一个列，标记需要显示标签的基因（例如top10最显著的）
ggplot(data = all_deg_results, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: AD vs HC",
       x = "log2 Fold Change (AD/HC)",
       y = "-log10(P-value)") +
  geom_vline(xintercept = c(-0, 0), linetype = "dashed") + # 添加logFC=0的线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") # 添加P=0.05的线

# # 火山图 (Volcano Plot) - 增强版
# all_deg_results$genelabels <- ""
# # 选择要标记的基因，例如按P值排序的前10个
# top_genes <- rownames(all_deg_results)[order(all_deg_results$P.Value)][1:10]
# all_deg_results$genelabels[rownames(all_deg_results) %in% top_genes] <- top_genes
# ggplot(all_deg_results, aes(x = logFC, y = -log10(P.Value))) +
#   geom_point(aes(color = ifelse(adj.P.Val < 0.05 & abs(logFC)>1, "Significant", "Not Significant")), alpha=0.6) +
#   scale_color_manual(values = c("grey", "red")) +
#   theme_minimal() +
#   geom_text_repel(aes(label = genelabels), # 使用ggrepel避免重叠
#                   box.padding = 0.5,
#                   point.padding = 0.1,
#                   force = 20,
#                   max.overlaps = Inf,
#                   size = 3) +
#   geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#   labs(title = "Volcano Plot",
#        x = "log2 Fold Change",
#        y = "-log10(P-value)",
#        color = "Status")
# 3.1 火山图 (Volcano Plot)------------------------------


# 3.2 热图 (Heatmap)-------------------------------------
# 获取显著DEGs的名称
deg_genes <- rownames(deg_significant)
# 从原始矩阵中提取这些基因的表达数据
heatmap_data <- expr_mat[deg_genes, ]
# 为热图添加样本分组注释
annotation_col <- data.frame(Group = Group_mat$Group)
rownames(annotation_col) <- colnames(heatmap_data)
# 绘制热图
pheatmap(heatmap_data[,1:dim(sample_metadata)[1]],
         scale = "row", # 按行（基因）进行Z-score标准化，使表达模式更清晰
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE, # 基因太多，不显示名字
         show_colnames = TRUE,
         annotation_col = annotation_col,
         color = colorRampPalette(c("blue", "white", "red"))(100), # 蓝-白-红配色
         main = "Heatmap of Significant DEGs (EXP vs HC)")
# 3.2 热图 (Heatmap)-------------------------------------


# 3.3 韦恩图 (Venn Diagram)------------------------------
# 获取交集基因进行后续分析
deg_list1 <- rownames(deg_significant[300:1800,]) # 例如 EXP vs HC 的全部DEGs
deg_list2 <- rownames(deg_significant[1250:2000,]) # 例如 Male vs Female 的DEGs
deg_list3 <- rownames(deg_significant[1601:2500,]) # 例如 Early AD vs Late AD 的DEGs
overlaps <- venn(list(ADvsHC = deg_list1, MvsF = deg_list2, EarlyvsLate = deg_list3), show.plot = FALSE)
# 获取特定交集的基因名
intersections <- attr(overlaps, "intersections")
# 查看所有交集的名称
names(intersections)
# 获取特定交集的基因
common_genes <- intersections$`ADvsHC:MvsF`
# 获取只在某一个列表中出现的基因
only_in_ad_hc <- setdiff(deg_list1, union(deg_list2, deg_list3))

# # 绘制三向韦恩图
# venn.plot <- venn.diagram(
#   x = list(
#     EXPvsHC = deg_list1,
#     MvsF = deg_list2,
#     EarlyvsLate = deg_list3
#   ),
#   filename = NULL, # 不直接保存到文件，先存为对象
#   category.names = c("AD vs HC", "Male vs Female", "Early vs Late"), # 设置类别名称
#   output = TRUE,
#   fill = c("cornflowerblue", "green", "yellow"), # 每个圈的填充色
#   alpha = 0.5, # 填充色透明度
#   cat.cex = 1.0, # 类别名称的字体大小
#   cat.dist = 0.05, # 类别名称到圈的距离
#   cat.pos = c(-30, 30, 180), # 类别名称的位置（角度）
#   margin = 0.1,
#   main = "Overlap of DEGs from Different Comparisons",
#   main.cex = 1.2
# );
# # 在R中显示图形
# try(dev.off(dev.list()["RStudioGD"]), silent=TRUE )
# grid.draw(venn.plot);
# # 保存图形为PNG
# png("Venn_Diagram_DEGs.png", width = 800, height = 800, res = 150)
# grid.draw(venn.plot)
# dev.off()
# # 3.3 韦恩图 (Venn Diagram)-------------------------------


# 3.4 UpSet图（对于比较多个集合）------------------------
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
# 3.4 UpSet图（对于比较多个集合）------------------------


# 3.5 MA图 (MA Plot)-------------------------------------
fit <- lmFit(expr_mat[,1:dim(sample_metadata)[1]], design_mat)
contrast_matrix <- makeContrasts(EXPvsHC = EXP - HC, levels = design_mat)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
plotMA(fit2, coef = "ADvsHC", main = "MA Plot: AD vs HC", status = NULL)
abline(h = 0, col = "red", lty = 2) # 添加logFC=0的参考线
# 可以自己着色显著基因
deg_status <- ifelse(all_deg_results$adj.P.Val < 0.05 & abs(all_deg_results$logFC) > 1, "Significant", "Not Sig")
plot(all_deg_results$AveExpr, all_deg_results$logFC,
     col = ifelse(deg_status == "Significant", "red", "black"),
     pch = 16, cex = 0.6,
     xlab = "Average Expression (A)", ylab = "log2 Fold Change (M)",
     main = "MA Plot")
abline(h = 0, col = "blue", lty = 2)
legend("topright", legend = c("Significant", "Not Significant"), col = c("red", "black"), pch = 16)
# 3.5 MA图 (MA Plot)-------------------------------------



### 4. 导出结果-----------------------------------------------------------------
# 导出所有基因的完整结果
file_All_Gene <- paste(output_dir, gse_number, '/All_Gene_Results_EXP_vs_HC.csv', sep = "", collapse = NULL)
write.csv(all_deg_results, file = file_All_Gene, row.names = TRUE)
# 导出显著的DEGs列表
file_DEGs <- paste(output_dir, gse_number, '/Significant_DEGs_EXP_vs_HC_Pval_0.05_logFC_0.csv', sep = "", collapse = NULL)
write.csv(deg_significant, file = file_DEGs, row.names = TRUE)
