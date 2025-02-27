# File: main.R 


# 初始化环境
rm(list = ls())
options(stringsAsFactors = FALSE)

# 设置路径
input_dir <- "raw_data"
output_dir <- "results"

# 加载子模块
source("src/diff_expr_analysis_script.R")
source("src/visualization_volcano_plot.R")

# 主分析流程
run_analysis <- function() {
  message("\n===== Starting Differential Expression Analysis =====")
  
  # 获取输入文件
  input_files <- list.files(input_dir, pattern = "matrix\\.txt$", full.names = TRUE)
  if (length(input_files) == 0) {
    stop("No input files found in: ", input_dir)
  }
  
  message("Found ", length(input_files), " input files")
  
  # 执行差异表达分析
  process_files(input_files)
  
  # 生成可视化结果
  message("\n===== Generating Volcano Plots =====")
  generate_all_volcano()
  
  message("\n===== Analysis Pipeline Completed Successfully =====")
}

# 执行主程序
system.time({
  tryCatch({
    run_analysis()
  }, error = function(e) {
    message("\n[ERROR] Pipeline failed: ", e$message)
    quit(status = 1)
  })
})
