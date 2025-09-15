# step2. Analysis DEGs
library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)

rm(list = ls()) #清空变量
gc()


### 1.数据读取----------------------------------------------
output_dir = './DATA/'
gse_number = 'GSE83500'
filename <-  paste(output_dir, gse_number, '/step1output.Rdata', sep = "", collapse = NULL)
load(file = filename)


### 2.准备数据结构---------------------------------------------
expr_mat <- combined_data
rownames(expr_mat) <- expr_mat$PROBEID
Group=ifelse(str_detect(sample_metadata$title,"non-MI"),
             "non-MI",
             "MI")
Group_ID <- paste0(Group, ave(Group, Group, FUN = seq_along))
Group_mat <- data.frame(sample_metadata$geo_accession,Group_ID = Group_ID,Group = Group)
colnames(Group_mat) <- c('GSM', 'SampleID', 'Group')
df2 <- data.frame(GSM=c('PROBEID','SYMBOL','ENTREZID'),
                  SampleID=c('PROBEID','SYMBOL','ENTREZID'),
                  Group=c('PROBEID','SYMBOL','ENTREZID'))
Group_mat <- rbind(Group_mat, df2)
colnames(expr_mat) <- Group_mat$SampleID


### 2.数据预处理---------------------------------------------
# 检查表达矩阵和样本信息是否匹配
stopifnot(colnames(expr_mat) == Group_mat$SampleID)
# 检查是否有缺失值或无限值
sum(is.na(expr_mat))
expr_mat_vector <- unlist(expr_mat)
# sum(!is.finite(x_vector)) # 应该都为0
sum(is.finite(expr_mat_vector))
# 可以绘制箱线图查看数据分布是否整齐（已归一化的数据应该很整齐）
boxplot(expr_mat[,1:dim(sample_metadata)[1]], main = "Sample Boxplots (after normalization)", las = 2)


### 3. 构建设计矩阵------------------------------------------
# 创建设计矩阵
# ~ 0 + Group 表示创建一个没有截距项的模型，直接为每个组生成一个系数
design_mat <- model.matrix(~ 0 + Group_mat$Group)
design_mat <- design_mat[1:dim(sample_metadata)[1], 2:3]
colnames(design_mat) <- c('EXP','HC') # 将列名命名为 "HC" 和 "AD"


### 4. 拟合线性模型------------------------------------------
# 将表达矩阵和设计矩阵拟合到线性模型中
fit <- lmFit(expr_mat[,1:dim(sample_metadata)[1]], design_mat)
# 创建对比矩阵，比较 EXP 相对于 HC 的变化
contrast_matrix <- makeContrasts(EXPvsHC = EXP - HC, levels = design_mat)
# 将对比矩阵应用到之前拟合的模型上
fit2 <- contrasts.fit(fit, contrast_matrix)
# 执行经验贝叶斯平滑，以获得更稳健的标准差估计
fit2 <- eBayes(fit2)
# 提取所有基因的统计结果
all_deg_results <- topTable(fit2, number = Inf, coef = "EXPvsHC", adjust.method = "BH", sort.by = "P")
# coef = "ADvsHC": 指定要提取的对比组
# adjust.method = "BH": 使用Benjamini & Hochberg方法计算FDR（调整后的p值）
# number = Inf: 提取所有基因
# sort.by = "P": 按P值排序

# 查看结果的前几行
head(all_deg_results)
# logFC: log2 Fold Change (AD/HC)。正数表示在AD中上调，负数表示下调。
# AveExpr: 平均表达量（所有样本的平均log2表达水平）。
# t: moderated t-statistic。
# P.Value: 原始的p值。
# adj.P.Val: 经过多重检验校正后的p值（FDR）。
# B: B统计量（log-odds）。




### 6.保存数据 ----------------------------------------------------
filename <-  paste(output_dir, gse_number, '/step2output.Rdata', sep = "", collapse = NULL)
save(gse_number, annot_table, sample_metadata, expr_mat, all_deg_results, design_mat, Group_mat, gpl_number, file = filename)
