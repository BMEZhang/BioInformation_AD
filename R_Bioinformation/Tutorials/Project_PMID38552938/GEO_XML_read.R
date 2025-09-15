# 安装并加载必要的包
if (!require("XML")) install.packages("XML")
if (!require("dplyr")) install.packages("dplyr") # 用于数据操作
if (!require("reshape2")) install.packages("reshape2") # 用于转换数据格式

library(XML)
library(dplyr)
library(reshape2)

# 设置文件路径
xml_file <- "path/to/your/GSE12345_family.xml" # 请替换为你的实际文件路径


# 解析XML文档
doc <- xmlParse(xml_file)




# 提取所有 'Sample' 节点
sample_nodes <- getNodeSet(doc, "//Sample")

# 创建一个空列表来存储每个样本的信息
sample_list <- list()

# 循环遍历每个 Sample 节点
for (i in 1:length(sample_nodes)) {
  sample <- sample_nodes[[i]]
  
  # 提取样本 accession (e.g., GSMXXXX)
  sample_id <- xmlAttrs(sample)["iid"]
  
  # 提取所有 'Channel' 下的 'Value' 标签内容
  # 这些通常是样本的特征，如 genotype, treatment, tissue 等
  characteristics <- xpathSApply(sample, ".//Value", xmlValue)
  # 获取 'Value' 标签的 'tag' 属性，作为列名
  characteristic_names <- xpathSApply(sample, ".//Value", xmlGetAttr, "tag")
  
  # 将特征值和名称组合成一个命名向量
  if (length(characteristics) > 0) {
    names(characteristics) <- characteristic_names
  }
  
  # 将样本ID和其他信息组合成一个数据框
  sample_df <- data.frame(sample_id = sample_id, 
                          as.list(characteristics), # 将命名向量转换为列表再转换为数据框的列
                          stringsAsFactors = FALSE)
  
  # 将这个样本的数据框添加到列表中
  sample_list[[i]] <- sample_df
}

# 将列表中的所有样本数据框合并成一个大的数据框
sample_metadata <- bind_rows(sample_list)

# 查看前几行
head(sample_metadata)




# 方法一：更手动的方法，理解结构
# 获取所有包含表达数据的 'Data-Table' 节点（通常每个样本一个）
data_tables <- getNodeSet(doc, "//Sample//Channel//Data-Table")

# 初始化一个列表来存储每个样本的表达向量
expression_list <- list()

for (i in 1:length(data_tables)) {
  sample_id <- xmlGetAttr(getNodeSet(doc, "//Sample")[[i]], "iid") # 获取对应样本的ID
  
  # 提取该 Data-Table 下所有 'Row' 的 'val' 属性（即表达值）
  expr_values <- as.numeric(xpathSApply(data_tables[[i]], ".//Row", xmlGetAttr, "val"))
  # 提取 'Row' 的 'idx' 属性（通常是探针ID），作为名字
  probe_ids <- xpathSApply(data_tables[[i]], ".//Row", xmlGetAttr, "idx")
  
  names(expr_values) <- probe_ids
  expression_list[[sample_id]] <- expr_values
}

# 将列表转换为矩阵，行是探针，列是样本
expression_matrix <- do.call(cbind, expression_list)
head(expression_matrix[, 1:3]) # 查看矩阵的前几行和前3列




# 提取样本元数据 (使用上面的代码)
# ... [插入提取 sample_metadata 的代码] ...

# 提取表达矩阵 (使用上面的代码)
# ... [插入提取 expression_matrix 的代码] ...

# 假设 sample_metadata 中有一列 'treatment' 包含分组信息（如 "Control", "Treated"）
# 假设我们想查看探针 "202783_at" 的表达情况

# 将矩阵转换为长格式，方便用ggplot2绘图
df_to_plot <- data.frame(
  sample_id = colnames(expression_matrix),
  expression = expression_matrix["202783_at", ], # 选择特定探针
  treatment = sample_metadata$treatment # 从元数据中添加分组信息
)

# 绘制盒形图
library(ggplot2)
ggplot(df_to_plot, aes(x = treatment, y = expression, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) + # 添加点显示样本分布
  labs(title = "Expression of Probe 202783_at",
       x = "Treatment Group",
       y = "Expression Level") +
  theme_minimal()

