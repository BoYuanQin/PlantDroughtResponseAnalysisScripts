# File: src/visualization_volcano_plot.R

library(dplyr)
library(ggplot2)
library(ggrepel)

generate_volcano_plot <- function(results_file, 
                                 output_dir = "results/visualization",
                                 top_n = 5) { 
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (!file.exists(results_file)) {
    stop("Results file not found: ", results_file)
  }
  
  data <- read.delim(results_file, stringsAsFactors = FALSE)
  sample_name <- tools::file_path_sans_ext(basename(results_file))
  
  # 准备绘图数据
  data <- data %>%
    mutate(
      logFDR = -log10(FDR),
      Significant = ifelse(FDR < 0.01 & abs(logFC) > 2, 
                          "Significant", "Not Significant")
    ) %>%
    # 添加标签逻辑
    arrange(FDR) %>% 
    mutate(
      Rank = row_number(),
      Label = ifelse(
        Significant == "Significant" & Rank <= top_n, 
        GeneID, 
        NA_character_
      )
    )
  
  # 创建火山图
   
  
  p <- ggplot(data, aes(x = logFC, y = logFDR, color = Significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("gray", "red")) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
    # 添加基因标签
    geom_text_repel(
      aes(label = Label),
      box.padding = 0.5,
      max.overlaps = Inf,
      segment.color = "grey50",
      size = 3,
      na.rm = TRUE  # 忽略NA值
    ) +
    labs(title = paste("Volcano Plot:", sample_name),
        x = "Log2 Fold Change",
        y = "-Log10 FDR") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "top",
      panel.grid.major = element_line(color = "grey90"),
      panel.border = element_rect(color = "black", fill = NA)
    )
  
  # 保存图形
  output_file <- file.path(output_dir, paste0(sample_name, "_volcano.pdf"))
  ggsave(output_file, p, width = 8, height = 6)
  message("Saved volcano plot with ", sum(!is.na(data$Label)), 
         " gene labels: ", output_file)
}


generate_all_volcano <- function(results_dir = "results/diff_expression") {
  dir.create("results/visualization", showWarnings = FALSE)
  
  result_files <- list.files(results_dir, 
                            pattern = "_full_results\\.tsv$", 
                            full.names = TRUE)
  
  if (length(result_files) == 0) {
    stop("No result files found in: ", results_dir)
  }
  
  lapply(result_files, function(f) {
    generate_volcano_plot(f)
  })
}
