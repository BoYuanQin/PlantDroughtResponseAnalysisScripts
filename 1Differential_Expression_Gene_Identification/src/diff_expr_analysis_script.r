# File: src/diff_expr_analysis_script.R 

process_files <- function(file_paths) {
  dir.create("results/diff_expression", showWarnings = FALSE, recursive = TRUE)
  
  library(parallel)

  cl <- makeCluster(detectCores() - 1)
  
  clusterEvalQ(cl, {
    library(methods)
    library(edgeR)
    library(readxl)
    
    
    HighSeqAnalyzer <- setRefClass(
      "HighSeqAnalyzer",
      fields = list(
        input_file_path = "character",
        significant_genes = "character",
        diff_expr_results = "data.frame"
      ),
      methods = list(
        initialize = function(input_file_path) {
          if (!file.exists(input_file_path)) {
            stop(paste("Input file not found:", input_file_path))
          }
          .self$input_file_path <- input_file_path
        },
        
        read_data = function(sep = "\t") {
          message("\nReading data from: ", .self$input_file_path)
          
          if (grepl("\\.txt$", .self$input_file_path)) {
            data <- tryCatch(
              read.table(.self$input_file_path, header = TRUE, sep = sep, stringsAsFactors = FALSE),
              error = function(e) stop("Error reading text file: ", e$message)
            )
          } else if (grepl("\\.xlsx$", .self$input_file_path)) {
            data <- tryCatch(
              read_excel(.self$input_file_path),
              error = function(e) stop("Error reading Excel file: ", e$message)
            )
          } else {
            stop("Unsupported file type: ", tools::file_ext(.self$input_file_path))
          }
          
          data <- as.data.frame(data)
          data <- na.omit(data)
          data <- aggregate(. ~ data[,1], data[,-1], mean)
          rownames(data) <- data[,1]
          data <- data[,-1, drop = FALSE]
          data[data < 0] <- 0
          as.matrix(data)
        },
        
        perform_differential_expression = function(data, sample_name) {
          message("Processing sample: ", sample_name)
          
          conditions <- ifelse(grepl("control", colnames(data), ignore.case = TRUE), 
                              "Control", "Treatment")
          
          # EdgeR分析流程
          y <- DGEList(counts = data, group = conditions)
          keep <- filterByExpr(y)
          y <- y[keep, , keep.lib.sizes = FALSE]
          y <- calcNormFactors(y, method = "TMM")
          
          # 保存过滤后的计数数据
          counts_data <- data.frame(GeneID = rownames(y$counts), y$counts)
          output_file <- file.path("results/diff_expression", 
                                  paste0(sample_name, "_filtered_counts.tsv"))
          write.table(counts_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
          
          # 差异表达分析
          design <- model.matrix(~0 + conditions)
          y <- estimateDisp(y, design)
          fit <- glmFit(y, design)
          lrt <- glmLRT(fit)
          diff_genes <- as.data.frame(topTags(lrt, n = Inf))
          
          # 保存完整结果
          full_results <- data.frame(GeneID = rownames(diff_genes), diff_genes)
          output_file <- file.path("results/diff_expression", 
                                  paste0(sample_name, "_full_results.tsv"))
          write.table(full_results, output_file, sep = "\t", row.names = FALSE)
          
          diff_genes
        },
        
        find_significant_genes = function(diff_matrix, sample_name) {
          sig_genes <- subset(diff_matrix, FDR < 0.01 & abs(logFC) > 2)
          sig_df <- data.frame(GeneID = rownames(sig_genes),
                              logFC = sig_genes$logFC,
                              FDR = sig_genes$FDR)
          
          output_file <- file.path("results/diff_expression",
                                  paste0(sample_name, "_significant_genes.tsv"))
          write.table(sig_df, output_file, sep = "\t", row.names = FALSE)
          
          rownames(sig_genes)
        },
        
        process = function(sample_name) {
          raw_data <- read_data()
          diff_results <- perform_differential_expression(raw_data, sample_name)
          .self$diff_expr_results <- diff_results
          .self$significant_genes <- find_significant_genes(diff_results, sample_name)
        }
      )
    )
  })

  clusterExport(cl, "file_paths", envir = environment()) 
  
  results <- parLapply(cl, seq_along(file_paths), function(i) {
    file_path <- file_paths[i]
    sample_name <- tools::file_path_sans_ext(basename(file_path))
    
    analyzer <- HighSeqAnalyzer$new(file_path)
    analyzer$process(sample_name)
    
    list(
      sample = sample_name,
      sig_genes = analyzer$significant_genes
    )
  })
  
  stopCluster(cl)
  
  # 提取共同差异基因
  all_genes <- lapply(results, function(x) x$sig_genes)
  common_genes <- Reduce(intersect, all_genes)
  writeLines(common_genes, 
            file.path("results/diff_expression", "common_significant_genes.txt"))
  
  message("\nAnalysis completed. Results saved to: results/diff_expression/")
  
  # 生成差异基因最多的样本矩阵
  select_top_sample <- function(diff_dir = "results/diff_expression") {
    sig_files <- list.files(diff_dir, pattern = "_significant_genes\\.tsv$", full.names = TRUE)
    
    if(length(sig_files) == 0) {
      stop("No significant genes files found in: ", diff_dir)
    }
    
    # 统计各样本差异基因数
    gene_counts <- sapply(sig_files, function(f) {
      nrow(read.delim(f))
    })
    
    # 获取最优样本
    top_sample <- names(which.max(gene_counts))
    top_sample <- tools::file_path_sans_ext(basename(top_sample))
    top_sample <- gsub("_significant_genes$", "", top_sample)
    
    message("Selected top sample: ", top_sample, " with ", max(gene_counts), " DEGs")
    return(top_sample)
  }
  
  # 生成合并计数矩阵
  generate_top_deg_matrix <- function(diff_dir = "results/diff_expression",
                                    output_dir = "results") {
    # 获取公共基因
    common_genes_file <- file.path(diff_dir, "common_significant_genes.txt")
    if(!file.exists(common_genes_file)) {
      stop("Common genes file not found: ", common_genes_file)
    }
    common_genes <- readLines(common_genes_file)
    
    # 自动选择最优样本
    top_sample <- select_top_sample(diff_dir)
    pattern <- paste0(top_sample, "_filtered_counts\\.tsv$")
    
    # 读取所有计数数据
    counts_files <- list.files(diff_dir, pattern = pattern, full.names = TRUE)
    if(length(counts_files) == 0) {
      stop("No counts files found with pattern: ", pattern)
    } 
    else if(length(counts_files) == 1) {
      # 单个文件直接读取
      deg_data_counts <- read.delim(counts_files[1], stringsAsFactors = FALSE)
    }
    
    # 筛选公共基因
    final_counts <- deg_data_counts %>%
      filter(GeneID %in% common_genes) %>%
      arrange(GeneID)
    
    # 保存结果
    output_file <- file.path(output_dir, "deg_data_counts.tsv")
    write.table(final_counts, output_file, sep = "\t", row.names = FALSE)
    message("\ndeg counts saved to: ", output_file)
  }  

  # 获取后续分析数据
  generate_top_deg_matrix()
  
  message("\nAnalysis completed. Results saved to: results/")
}
