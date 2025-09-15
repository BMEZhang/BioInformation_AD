# step 1
library(GEOquery)
library(Biobase)
library(stringr)
library(dplyr)
library(AnnotationDbi)
library(limma)

rm(list = ls()) #清空变量
gc()


###  1.数据读取----------------------------------------------
output_dir = './DATA/'
gse_number = 'GSE83500'
filename <-  paste(output_dir, gse_number, '/step0output.Rdata', sep = "", collapse = NULL)
load(file = filename)
# 1.Group----实验分组要去阅读临床信息的表格来获取，每个GSE的都不一样
# Group=ifelse(str_detect(sample_metadata$source_name_ch1,"non-MI"),
#              "non-MI",
#              "MI")
#设置参考水平，指定levels，对照组在前，处理组在后
# Group = factor(Group,
#                levels = c("non-MI","MI"))
# 只想改变参考水平而不改变其他水平的顺序
# Group <- relevel(Group, ref = "non-MI")


### 2.探针注释的获取-----------------------------------------
# 2.1 获取探针注释包：优先使用平台专用注释包
# 2.2 处理多对一映射：合理处理多个探针对应同一基因的情况
# 2.3 验证注释质量：检查注释覆盖率和唯一基因数量
# 2.4 保持一致性：在整个分析流程中使用相同的注释方法和标准
# 2.5 记录处理步骤：确保分析的可重复性

# 2.1 获取探针注释包
gpl <- getGEO(gpl_number, destdir = ".")  # 下载平台信息
annot_table <- Table(gpl)   # 提取注释表
colnames(annot_table)       # 查看可用列
probe_annot <- annot_table[, c("ID", "Gene Symbol", "Entrez Gene")]  # 选择需要的列（通常是ID和Gene Symbol）列名可能因平台而异
colnames(probe_annot) <- c("PROBEID", "SYMBOL", "ENTREZID")  # 标准化列名
# 过滤掉没有注释到基因的探针（NA值）
probe_annot <- probe_annot[!is.na(probe_annot$SYMBOL) & !is.na(probe_annot$ENTREZID), ] 

# 2.2 验证注释质量
# 检查注释覆盖率
probe_ids <- rownames(probe_annot)
total_probes <- length(probe_ids)
annotated_probes <- sum(!is.na(probe_annot$SYMBOL) & probe_annot$SYMBOL != "")
coverage_rate <- annotated_probes / total_probes * 100
cat(sprintf("注释覆盖率: %.2f%%\n", coverage_rate))
# 检查基因数量
unique_genes <- length(unique(probe_annot$SYMBOL[!is.na(probe_annot$SYMBOL)]))
cat(sprintf("唯一基因数量: %d\n", unique_genes))

# # 2.2 处理多对一映射 取平均值
# probe_ids <- rownames(sample_matrix)
# # 计算平均表达量并选择最高表达的探针
# mean_expr <- rowMeans(sample_matrix)
# probe_annot$mean_expr <- mean_expr[probe_annot$PROBEID]
# # 确保所有列为普通向量类型
# probe_annot <- as.data.frame(lapply(probe_annot, as.vector), stringsAsFactors = FALSE)
# # 按基因分组并找到每个基因中表达量最高的探针
# probe2gene_ordered <- probe_annot[order(probe_annot$SYMBOL, -probe_annot$mean_expr), ]
# # 使用 duplicated() 函数获取每个基因的第一个（最高表达）探针
# selected_probes <- probe2gene_ordered[!duplicated(probe2gene_ordered$SYMBOL), "PROBEID"]
# # 筛选表达矩阵
# expr_filtered <- sample_matrix[selected_probes, ]
# rownames(expr_filtered) <- probe2gene_ordered$SYMBOL[match(selected_probes, probe2gene_ordered$PROBEID)]


### 3 舍弃多对一映射探针-------------------------------------------
# 3.1 找出有多映射的探针
multi_mapping_probes <- probe_annot$PROBEID[grepl("///", probe_annot$ENTREZID)]
# 3.2 从表达矩阵中过滤掉这些探针
expr_mat_filtered <- sample_matrix[!rownames(sample_matrix) %in% multi_mapping_probes, ]
# 3.3 从注释矩阵中过滤掉这些探针
probe_annot_filtered <- probe_annot[!probe_annot$PROBEID %in% multi_mapping_probes, ]
# 3.4 从注释矩阵中过滤掉这些探针
# 确保表达矩阵是一个数据框（Dataframe）
expr_df <- as.data.frame(expr_mat_filtered)
expr_df$PROBEID <- rownames(expr_df) # 将行名变为一列
# 3.5 使用dplyr的left_join根据PROBEID列合并表达矩阵和注释信息
combined_data <- left_join(expr_df, probe_annot_filtered, by = "PROBEID")
probe_ids <- rownames(combined_data) 

### 4.保存数据 ----------------------------------------------------
filename <-  paste(output_dir, gse_number, '/step1output.Rdata', sep = "", collapse = NULL)
save(gse_number,sample_metadata,annot_table,combined_data,gpl_number,file = filename)



# # 2.3 创建注释感知的表达矩阵
# annotated_expr_matrix <- sample_matrix %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column("PROBEID") %>%
#   left_join(annotations, by = "PROBEID") %>%
#   # 去除没有注释的探针
#   filter(!is.na(SYMBOL) & SYMBOL != "") %>%
#   # 可以选择保留一个探针 per 基因
#   group_by(SYMBOL) %>%
#   slice(1) %>%  # 保留每个基因的第一个探针
#   ungroup() %>%
#   # 将基因名作为行名
#   tibble::column_to_rownames("SYMBOL") %>%
#   select(-PROBEID, -ENTREZID)  # 移除不需要的列
# 
# # 转换为数值矩阵
# annotated_expr_matrix <- as.matrix(annotated_expr_matrix)




