# File: src/visualization_volcano_plot.R

library(dplyr)
library(ggplot2)
library(ggrepel)

generate_volcano_plot <- function(merged_results, common_genes,
                                 output_dir = "results/visualization",
                                 top_n = 10) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 标记显著性并区分上下调
  plot_data <- merged_results %>%
    mutate(
      Significance = case_when(
        GeneID %in% common_genes & logFC > 0 ~ "Upregulated",
        GeneID %in% common_genes & logFC < 0 ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      logFDR = -log10(FDR)
    ) %>%
    arrange(FDR) %>%
    group_by(Significance) %>%
    mutate(
      Rank = ifelse(Significance != "Not Significant", row_number(), NA_integer_),
      Label = ifelse(Rank <= top_n, GeneID, NA_character_)
    ) %>%
    ungroup()
  
  # 自定义颜色方案
  color_values <- c("Upregulated" = "#E64B35", 
                   "Downregulated" = "#3182BD", 
                   "Not Significant" = "grey60")
  
  p <- ggplot(plot_data, aes(x = logFC, y = logFDR, color = Significance)) +
    geom_point(alpha = 0.6, size = 1.8) +
    scale_color_manual(values = color_values) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey30") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    geom_text_repel(
      aes(label = Label),
      box.padding = 0.6,
      max.overlaps = 20,
      segment.color = "grey50",
      size = 3.2,
      na.rm = TRUE
    ) +
    labs(
      title = "Combined Volcano Plot (Averaged Results)",
      x = expression(log[2]("Fold Change")),
      y = expression(-log[10]("FDR"))
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "grey50"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2)
    )
  
  output_file <- file.path(output_dir, "combined_volcano_plot.pdf")
  ggsave(output_file, p, width = 10, height = 7)
  message("Saved combined volcano plot to: ", output_file)
}

generate_all_volcano <- function(results_dir = "results/diff_expression") {
  # 读取所有样本结果文件并合并
  result_files <- list.files(results_dir, 
                            pattern = "_full_results\\.tsv$", 
                            full.names = TRUE)
  
  if (length(result_files) == 0) {
    stop("No result files found in: ", results_dir)
  }
  
  # 合并数据并取均值
  merged_results <- lapply(result_files, function(f) {
    read.delim(f) %>% select(GeneID, logFC, FDR)
  }) %>% 
    bind_rows() %>%
    group_by(GeneID) %>%
    summarise(
      logFC = mean(logFC, na.rm = TRUE),
      FDR = mean(FDR, na.rm = TRUE)
    )
  merged_results_path <- file.path(results_dir, "merged_results.tsv")
  write.table(merged_results, 
             file = merged_results_path,
             sep = "\t",          
             row.names = FALSE,   
             col.names = TRUE,    
             quote = FALSE,       
             na = "NA",           
             dec = ".")           
  message("Saved merged results to: ", merged_results_path)
  
  # 读取公共差异基因
  common_genes <- readLines(file.path(results_dir, "common_significant_genes.txt"))
  
  # 生成火山图
  generate_volcano_plot(merged_results, common_genes)
}
