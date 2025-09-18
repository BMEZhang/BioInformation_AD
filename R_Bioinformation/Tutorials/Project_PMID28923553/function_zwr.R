### 调用的函数



### library---------------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(GEOquery)
library(tidyr)
library(survival)
library(Biobase)
library(stringr)
library(DESeq2)
library(AnnotationDbi)
library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(gplots)
library(grid)
library(VennDiagram)
library(UpSetR)
library(ggrepel) # 防止标签重叠
# ------------------------------------------------------------------------------



### function--------------------------------------------------------------------
# 数据读取
data_download_zwr <- function(output_dir,gse_number){
  if (!file.exists(output_dir)) {dir.create(output_dir)}
  eSet <- getGEO(gse_number, 
                 destdir = output_dir, 
                 AnnotGPL = F,
                 getGPL = F)
  sample_matrix <- exprs(eSet[[1]]) # 提取表达矩阵
  sample_matrix <- examine_matrix_zwr(sample_matrix) # 检查数据是否已经标准化
  sample_metadata <- pData(eSet[[1]]) # 提取样本元数据
  gpl_number <- eSet[[1]]@annotation
  gpl <- getGEO(gpl_number,destdir = output_dir)
  annot_table <- Table(gpl)   # 提取原始注释表
  
  p = identical(rownames(sample_metadata),colnames(sample_matrix));p
  if(p){return(list(sample_matrix,sample_metadata,annot_table))}else{print('DOWNLOAD ERROR')}
  
}

# 提取PROBEID、SYMBOL以及ENTREZID
probe_annot_zwr <- function(annot_table,TABLE_ID,TABLE_SYMBOL,TABLE_ENTREZID){
  if ((TABLE_ENTREZID %in% colnames(annot_table)) & (TABLE_SYMBOL %in% colnames(annot_table))) {
    print("使用GPL平台注释文件")
    probe2entrez <- annot_table %>%
      dplyr::select(PROBEID = TABLE_ID, SYMBOL = TABLE_SYMBOL, ENTREZID = TABLE_ENTREZID) %>%
      filter(!is.na(ENTREZID) & ENTREZID != "")
    
    multi_mapping_probes <- probe2entrez$SYMBOL[grepl("///", probe2entrez$ENTREZID)]
    probe_annot_filtered <- probe2entrez[!probe2entrez$SYMBOL %in% multi_mapping_probes, ]
    probe_annot_filtered <- probe_annot_filtered %>%
      mutate(
        ENTREZID = na_if(as.character(ENTREZID), ""),     # 将空字符串 "" 转换为 NA
        ENTREZID = na_if(as.character(ENTREZID), "---"),  # 将 "---" 转换为 NA
        ENTREZID = na_if(as.character(ENTREZID), "unknown") # 如果需要，还可以处理其他标记
      ) %>%
      na.omit() # 现在，na.omit() 会一次性删除所有包含NA的行
    
    multi_mapping_probes <- probe2entrez$ENTREZID[grepl("///", probe2entrez$SYMBOL)]
    probe_annot_filtered <- probe2entrez[!probe2entrez$ENTREZID %in% multi_mapping_probes, ]
    probe_annot_filtered <- probe_annot_filtered %>%
      mutate(
        SYMBOL = na_if(SYMBOL, ""),     # 将空字符串 "" 转换为 NA
        SYMBOL = na_if(SYMBOL, "---"),  # 将 "---" 转换为 NA
        SYMBOL = na_if(SYMBOL, "unknown") # 如果需要，还可以处理其他标记
      ) %>%
      na.omit() # 现在，na.omit() 会一次性删除所有包含NA的行
    
    examine_gpl_zwr(probe_annot_filtered) # 检查GPL去空后数据
    dup <- duplicated(probe_annot_filtered[, 2])
    probe_annot <- probe_annot_filtered[!dup, ]
    # probe_annot <- probe_annot_filtered[unique(probe_annot_filtered$SYMBOL),]
  }else if (TABLE_SYMBOL %in% colnames(annot_table))
  {
    print("缺少ENTREZID，使用org.Hs.eg.db注释文件")
    probe2entrez <- annot_table %>%
      dplyr::select(PROBEID = TABLE_ID, SYMBOL = TABLE_SYMBOL) %>%
      filter(!is.na(PROBEID) & SYMBOL != "")
    probe2entrez <- probe2entrez[!duplicated(probe2entrez[, "SYMBOL"]), ]
    SYMBOL_ids <- probe2entrez$SYMBOL
    id_mapping <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = SYMBOL_ids,    # 你的探针ID向量
                                        keytype = "SYMBOL",  # 
                                        columns = c("SYMBOL", "ENTREZID"))
    id_mapping <- id_mapping[!duplicated(id_mapping[, "SYMBOL"]), ]
    probe_annot_filtered <- merge(probe2entrez, id_mapping, by = "SYMBOL")
    probe_annot_filtered <- probe_annot_filtered %>%
      mutate(
        ENTREZID = na_if(ENTREZID, ""),     # 将空字符串 "" 转换为 NA
        ENTREZID = na_if(ENTREZID, "---"),  # 将 "---" 转换为 NA
        ENTREZID = na_if(ENTREZID, "unknown") # 如果需要，还可以处理其他标记
      ) %>%
      na.omit() # 现在，na.omit() 会一次性删除所有包含NA的行
    examine_gpl_zwr(probe_annot_filtered) # 检查GPL去空后数据
    dup <- duplicated(probe_annot_filtered[, "SYMBOL"])
    probe_annot <- probe_annot_filtered[!dup, ]
  }else
  {
    print("缺少SYMBOL，使用org.Hs.eg.db注释文件")
    probe2entrez <- annot_table %>%
      dplyr::select(PROBEID = TABLE_ID, ENTREZID = TABLE_ENTREZID) %>%
      filter(!is.na(PROBEID) & ENTREZID != "")
    probe2entrez <- probe2entrez[!duplicated(probe2entrez[, "ENTREZID"]), ]
    ENTREZID_ids <- probe2entrez$ENTREZID
    id_mapping <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = as.character(ENTREZID_ids),    # 你的探针ID向量
                                        keytype = "ENTREZID",  # 
                                        columns = c("SYMBOL", "ENTREZID"))
    id_mapping <- id_mapping[!duplicated(id_mapping[, "ENTREZID"]), ]
    probe_annot_filtered <- merge(probe2entrez, id_mapping, by = "ENTREZID")
    probe_annot_filtered <- probe_annot_filtered %>%
      mutate(
        SYMBOL = na_if(SYMBOL, ""),     # 将空字符串 "" 转换为 NA
        SYMBOL = na_if(SYMBOL, "---"),  # 将 "---" 转换为 NA
        SYMBOL = na_if(SYMBOL, "unknown") # 如果需要，还可以处理其他标记
      ) %>%
      na.omit() # 现在，na.omit() 会一次性删除所有包含NA的行
    examine_gpl_zwr(probe_annot_filtered) # 检查GPL去空后数据
    dup <- duplicated(probe_annot_filtered[, "ENTREZID"])
    probe_annot <- probe_annot_filtered[!dup, ]
    }
  return(probe_annot)
}

# 检查数据是否已经标准化
examine_matrix_zwr <- function(exp) {
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
# 检查GPL质量
examine_gpl_zwr <- function(probe_annot){
  probe_ids <- probe_annot$PROBEID
  total_probes <- length(probe_ids)
  annotated_symbol <- sum(!is.na(probe_annot$SYMBOL) & probe_annot$SYMBOL != "")
  coverage_rate <- annotated_symbol / total_probes * 100
  cat(sprintf("注释覆盖率: %.4f%%\n", coverage_rate))
  # 检查基因数量
  unique_symbol <- length(unique(probe_annot$SYMBOL[!is.na(probe_annot$SYMBOL)]))
  cat(sprintf("唯一SYMBOL数量: %d\n", unique_symbol))
  unique_genes <- length(unique(probe_annot$ENTREZID[!is.na(probe_annot$ENTREZID)]))
  cat(sprintf("唯一ENTREZID数量: %d\n", unique_genes))
  }

# 检查表达矩阵和样本信息是否匹配
examine_meta_zwr <- function(expr_mat, Group_mat){
  stopifnot(colnames(expr_mat) == Group_mat$SampleID)
  # 检查是否有缺失值或无限值
  temp <- sum(is.na(expr_mat))
  if (!temp) {print('expr_mat 没有缺失值')} else {print('expr_mat有缺失值,数量为: %d', temp)}
  expr_mat_vector <- unlist(expr_mat)
  # sum(!is.finite(x_vector)) # 应该都为0
  temp <- sum(!is.finite(expr_mat_vector))
  if (!temp) {print('expr_mat 没有极大值')} else {print('expr_mat有极大值,数量为: %d', temp)}
  # 可以绘制箱线图查看数据分布是否整齐（已归一化的数据应该很整齐）
  boxplot(expr_mat[,1:dim(expr_mat)[2]], main = "Sample Boxplots (after normalization)", las = 2)
  print('已绘制boxplot图')
}
# 舍弃一对多以及多对一的SYMBOL
desert_expr_zwr <- function(sample_matrix, probe_annot){
  # 处理可能存在的NA值（没有对应基因的探针）
  probe_annot <- probe_annot %>%
    mutate(
      SYMBOL = na_if(SYMBOL, ""),     # 将空字符串 "" 转换为 NA
      SYMBOL = na_if(SYMBOL, "---"),  # 将 "---" 转换为 NA
      SYMBOL = na_if(SYMBOL, "unknown") # 如果需要，还可以处理其他标记
    ) %>%
    na.omit() # 现在，na.omit() 会一次性删除所有包含NA的行
  # 1 找出有多映射的探针
  multi_mapping_probes <- probe_annot$PROBEID[grepl("///", probe_annot$ENTREZID)]
  # 2 从表达矩阵中过滤掉这些探针
  expr_mat_filtered <- sample_matrix[!rownames(sample_matrix) %in% multi_mapping_probes, ]
  # 3 从注释矩阵中过滤掉这些探针
  probe_annot_filtered <- probe_annot[!probe_annot$PROBEID %in% multi_mapping_probes, ]
  
  # 4 找出有多映射的基因
  # 4.1 确保表达矩阵是一个数据框（Dataframe），并将探针ID作为一列
  expr_df <- as.data.frame(expr_mat_filtered)
  expr_df$probe_id <- rownames(expr_df)
  # 4.2 将表达矩阵与探针注释信息合并
  merged_data <- merge(probe_annot_filtered, expr_df, by.x = "PROBEID", by.y = "probe_id")
  # 4.3. 分组并求平均值
  # 按基因符号(SYMBOL)分组，并对所有样本列(这里假设样本名列如 GSM1, GSM2...)求平均值
  gene_expr <- merged_data %>%
    dplyr::select(-PROBEID,-ENTREZID) %>%        # 移除探针ID列
    group_by(SYMBOL) %>%               # 按基因符号分组
    summarise(across(everything(), function(x) mean(x, na.rm = TRUE))) # 对所有其他列（样本）求平均值
  gene_expr <- as.data.frame(gene_expr)
  rownames(gene_expr) <- gene_expr$SYMBOL
  # row.names(gene_expr) <- gene_expr$SYMBOL
  # expr_df$PROBEID <- rownames(expr_df) # 将行名变为一列
  # 3.5 使用dplyr的left_join根据PROBEID列合并表达矩阵和注释信息
  # combined_data <- left_join(expr_df, probe_annot_filtered, by = "PROBEID")
  gene_expr <- select(gene_expr, -SYMBOL) # 删除 SYMBOL 列
  return(gene_expr)
}

# 保存expr_mat，sample_metadata，probe_annot至csv表以及.Rdata文件。
save_mat_zwr <- function(output_dir1, gse_number,expr_mat, sample_metadata,probe_annot){
  new_line <- rownames(expr_mat)
  expr_mat_file <- data.frame(expr_mat,SYMBOL = new_line)
  expr_mat_file = select(expr_mat_file,dim(expr_mat_file)[2],1:dim(expr_mat_file)[2]-1)
  exprdata_Name <- paste(output_dir1, gse_number, "_matrix_expr.csv", sep = "", collapse = NULL)
  write.csv(expr_mat_file, exprdata_Name, row.names = F)
  rm(expr_mat_file,new_line)
  cat(sprintf("%s文件已保存\n", exprdata_Name))
  
  new_line <- rownames(sample_metadata)
  sample_metadata_file <- data.frame(sample_metadata,SAMPLE = new_line)
  sample_metadata_file = select(sample_metadata_file,dim(sample_metadata_file)[2],1:dim(sample_metadata_file)[2]-1)
  metadata_Name <- paste(output_dir1, gse_number, "_metadata.csv", sep = "", collapse = NULL)
  write.csv(sample_metadata_file, metadata_Name, row.names = F)
  rm(sample_metadata_file,new_line)
  cat(sprintf("%s文件已保存\n", metadata_Name))
  
  filename <-  paste(output_dir1, gse_number, '_step1output.Rdata', sep = "", collapse = NULL)
  gpl_number <- unique(sample_metadata$platform_id)
  save(gse_number,probe_annot,sample_metadata,annot_table,expr_mat,gpl_number,file = filename) # 保存数据
  cat(sprintf("%s文件已保存\n", filename))
  
  png(file = paste(output_dir1, gse_number, '_Expr_boxplot.png', sep = "", collapse = NULL), width = 800, height = 600)
  boxplot(expr_mat, main = gse_number, ylab="Express",xlab="Sample")
  dev.off()
}


# 模型线性拟合
fit_line_zwr <- function(Group_mat,expr_mat){
  # 创建设计矩阵
  # ~ 0 + Group 表示创建一个没有截距项的模型，直接为每个组生成一个系数
  design_mat <- model.matrix(~ 0 + Group_mat$Group)
  design_mat <- design_mat[1:dim(expr_mat)[2], 1:2]
  colnames(design_mat) <- c('EXP','HC') # 将列名命名为 "EXP" 和 "HC"
  # 将表达矩阵和设计矩阵拟合到线性模型中
  fit <- lmFit(expr_mat[,1:dim(expr_mat)[2]], design_mat)
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
  return(all_deg_results)
}
# 结果可视化
vis_GOKEGG_zwr <- function(val,ego_bp,kk,output_dir1,gse_number){
  switch(EXPR = val,
         "barplot"= {
           png(file = paste(output_dir1, gse_number, '_GO_barplot.png', sep = "", collapse = NULL), width = 800, height = 600)
           barplot(ego_bp,
                   showCategory = 15, # 显示最显著的前15个Term
                   title = "GO Biological Process Enrichment",
                   font.size = 8)
           dev.off()
         },
         "dotplot"= {
           png(file = paste(output_dir1, gse_number, '_KEGG_dotplot.png', sep = "", collapse = NULL), width = 800, height = 600)
           dotplot(ego_bp, showCategory = 15)
           dotplot(kk,
                   showCategory = 15,
                   title = "KEGG Pathway Enrichment")
           dev.off()
         },
         "ego_bp2"= {
           png(file = paste(output_dir1, gse_number, '_GO_bp2.png', sep = "", collapse = NULL), width = 800, height = 600)
           ego_bp2 <- pairwise_termsim(ego_bp)
           emapplot(ego_bp2, showCategory = 20)
           dev.off()
         },
         "cnetplot"= {
           png(file = paste(output_dir1, gse_number, '_GO_cnetplot.png', sep = "", collapse = NULL), width = 800, height = 600)
           cnetplot(ego_bp,
                    showCategory = 5, # 只显示前5个Term，否则图会很乱
                    circular = FALSE,
                    colorEdge = TRUE,
                    node_label = "all")
           dev.off()
         }, # 显示所有节点标签
         {
           stop("Unknown plot type: ", val)  # 未知类型报错
         }
  )
}
# PPI分析
PPI_analysis_zwr <- function(converted){
  library(STRINGdb)
  # 创建一个STRINGdb对象（针对人类物种）
  string_db <- STRINGdb$new(version = "11.5", 
                            species = 9606, # 人的物种编号
                            score_threshold = 400, # 互作置信度阈值(0-1000)
                            input_directory = "")
  
  # 映射我们的基因（使用之前转换好的ENTREZID数据框）
  mapped_genes <- string_db$map(converted, "ENTREZID", removeUnmappedRows = TRUE)
  # mapped_genes <- string_db$map(converted, "ENSEMBL", removeUnmappedRows = TRUE)
  
  # 获取互作对（边列表）
  ppi_edges <- string_db$get_interactions(mapped_genes$STRING_id)
  
  # 为边列表添加原始的Gene Symbol信息
  # 首先创建一个映射字典
  string_to_symbol <- setNames(mapped_genes$SYMBOL, mapped_genes$STRING_id)
  ppi_edges$from_symbol <- string_to_symbol[ppi_edges$from]
  ppi_edges$to_symbol <- string_to_symbol[ppi_edges$to]
  return(ppi_edges)
  # 保存PPI边列表，用于导入Cytoscape或其他软件
  # write.csv(ppi_edges, "PPI_Edge_List_for_Cytoscape.csv", row.names = FALSE)
}
# 火山图绘制
Volcano_zwr <- function(all_deg_results,P_value=0.05,logFC_value=0){
  # 为绘图数据添加一列，标记显著基因
  all_deg_results$Significant <- ifelse(all_deg_results$P.Value < P_value & abs(all_deg_results$logFC) > logFC_value, "Significant", "Not Sig")
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
}
# 热图绘制
HeatMap_zwr <- function(deg_significant,expr_mat,Group_mat){
  # 获取显著DEGs的名称
  deg_genes <- rownames(deg_significant)
  # 从原始矩阵中提取这些基因的表达数据
  heatmap_data <- expr_mat[deg_genes, ]
  # 为热图添加样本分组注释
  annotation_col <- data.frame(Group = Group_mat$Group)
  rownames(annotation_col) <- colnames(heatmap_data)
  # 绘制热图
  pheatmap(heatmap_data[,1:dim(expr_mat)[2]],
           scale = "row", # 按行（基因）进行Z-score标准化，使表达模式更清晰
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = FALSE, # 基因太多，不显示名字
           show_colnames = TRUE,
           annotation_col = annotation_col,
           color = colorRampPalette(c("blue", "white", "red"))(100), # 蓝-白-红配色
           main = "Heatmap of Significant DEGs (EXP vs HC)")
  print('绘制heatmap')
}
# Ma图绘制
Ma_zwr <- function(Group_mat, expr_mat){
  # 创建设计矩阵
  # ~ 0 + Group 表示创建一个没有截距项的模型，直接为每个组生成一个系数
  design_mat <- model.matrix(~ 0 + Group_mat$Group)
  design_mat <- design_mat[1:dim(expr_mat)[2], 1:2]
  colnames(design_mat) <- c('EXP','HC') # 将列名命名为 "EXP" 和 "HC"
  fit <- lmFit(expr_mat[,1:dim(expr_mat)[2]], design_mat)
  contrast_matrix <- makeContrasts(EXPvsHC = EXP - HC, levels = design_mat)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  plotMA(fit2, coef = "EXPvsHC", main = "MA Plot: EXP vs HC", status = NULL)
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
}

# 韦恩图绘制
Venn_zwr <- function(deg_significant){
  # 获取数据
  deg_list1 <- deg_significant[[1]]
  deg_list2 <- deg_significant[[2]]
  deg_list3 <- deg_significant[[3]]
  # 绘制三向韦恩图
  venn.plot <- venn.diagram(
    x = list(
      EXPvsHC = deg_list1,
      MvsF = deg_list2,
      EarlyvsLate = deg_list3
    ),
    filename = NULL, # 不直接保存到文件，先存为对象
    category.names = c("AD vs HC", "Male vs Female", "Early vs Late"), # 设置类别名称
    output = TRUE,
    fill = c("cornflowerblue", "green", "yellow"), # 每个圈的填充色
    alpha = 0.5, # 填充色透明度
    cat.cex = 1.0, # 类别名称的字体大小
    cat.dist = 0.05, # 类别名称到圈的距离
    cat.pos = c(-30, 30, 180), # 类别名称的位置（角度）
    margin = 0.1,
    main = "Overlap of DEGs from Different Comparisons",
    main.cex = 1.2
  );
  # 在R中显示图形
  try(dev.off(dev.list()["RStudioGD"]), silent=TRUE )
  grid.draw(venn.plot);
  # 保存图形为PNG
  png("Venn_Diagram_DEGs.png", width = 800, height = 800, res = 150)
  grid.draw(venn.plot)
  dev.off()
}
# ------------------------------------------------------------------------------
