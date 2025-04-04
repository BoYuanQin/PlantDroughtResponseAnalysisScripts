# Rscript Rscript/Differential_Gene_Expression_Analysis.R GSE202931_matrix.txt GSE201380_matrix.txt GSE201379_matrix.txt 
HightSeqAnalyzer <- setRefClass(
  "HightSeqAnalyzer",
    fields = list(
      input_file_path  = "character",
      significant_genes = "character",
      diff_expr_results = "data.frame"
    ),
    methods = list(
      initialize = function(input_file_path) {
        .self$input_file_path <- input_file_path
      },
      read_data = function(file_path, sep = "\t") {
        library(readxl)
        if (grepl("\\.txt$", file_path)) {
          data <- read.table(file_path, header = TRUE, sep = sep,stringsAsFactors=FALSE)
        } else if (grepl("\\.xlsx$", file_path)) {
          data <- read_excel(file_path)
        } else {
          stop("Unsupported file type")
        }
        data <- as.data.frame(data)
        data <- na.omit(data)
        data <- aggregate(. ~ data[,1], data[,-1], mean)
        rownames(data) <- data[,1]
        data <- data[,-1]
        data[data < 0] <- 0
        data <- as.matrix(data)      
        return(data)
      },
      perform_differential_expression = function(data,file_name,output_dir) {
        library(edgeR)
        colnames_data <- colnames(data)
        conditions <- ifelse(grepl("control", colnames_data), "Control", "Treatment")
        y <- DGEList(counts = data, group = conditions)  
        keep <- filterByExpr(y, group = as.factor(conditions))
        y <- y[keep, , keep.lib.sizes = FALSE]
		counts_data <- y$counts
		# counts_data <- log1p(counts_data)
		counts_data <- data.frame(GeneId = rownames(counts_data), counts_data)
		output_file_path <- file.path(output_dir, paste(file_name, "AT_Keep_Genes_Expression.txt", sep="_"))
		write.table(counts_data, output_file_path, sep = "\t", row.names = FALSE,quote = FALSE)
		y <- calcNormFactors(y, method = "TMM")
        design <- model.matrix(~0 + conditions)
        colnames(design) <- levels(conditions)
		# print(head(y$counts, 10))
		y <- estimateDisp(y, design)
		fit <- glmFit(y, design)
		lrt <- glmLRT(fit, coef = 2)
		diff_genes <- as.data.frame(topTags(lrt, n=Inf))
		print(head(diff_genes,10))
        return(diff_genes)
      },
      find_significant_genes = function(diff_matrix,file_name,output_dir) {
        selected_genes <- subset(diff_matrix, PValue < 0.01 & abs(logFC) > 2)
		output_file_path <- file.path(output_dir, file_name, "AT_significant_genes.txt")
		selected_genes_expression <- data.frame(GeneId=rownames(selected_genes),selected_genes)
		write.table(selected_genes_expression,output_file_path, sep = "\t", row.names = FALSE)
        return(as.character(rownames(selected_genes)))
      },
      process = function(file_name,output_dir) {
        data <- read_data(.self$input_file_path)
        .self$diff_expr_results <- perform_differential_expression(data,file_name,output_dir)
        .self$significant_genes <- find_significant_genes(.self$diff_expr_results,file_name,output_dir)
      }
    )
)
#基因差异表达分析
process_files <- function(files, output_dir,input_dir) {
  library(parallel)
  no_cores <- detectCores() - 1
  results <- mclapply(files, function(file) {
    input_file_path <- file.path(input_dir, file)
	output_file_path <- file.path(output_dir, gsub("_matrix\\.txt$", "", basename(file)))
	file_name <- gsub("_matrix\\.txt$","",basename(file))
    dir.create(output_file_path, recursive = TRUE)
    analyzer <- HightSeqAnalyzer$new(input_file_path)
    analyzer$process(file_name,output_dir)
    return(list(file = file, significant_genes = analyzer$significant_genes))
  }, mc.cores = no_cores)
  all_significant_genes <- lapply(results, function(res) res$significant_genes)
  common_genes <- Reduce(intersect, all_significant_genes)
  write.table(common_genes, file.path(output_dir, "AT_common_significant_genes.txt"), sep = "\t", row.names = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Insufficient arguments provided. Usage: Rscript script.R <file1> [<file2> ...]", call. = FALSE)
}
files <- args
output_dir <- "Results"
input_dir <- "Data/拟南芥"
process_files(files, output_dir,input_dir)

#生成后期需要用到的基因表达矩阵
read_counts_data <- function(directory, pattern) {
  file_paths <- list.files(directory, pattern = pattern, full.names = TRUE)
  list_of_counts_data <- lapply(file_paths, function(file_path) {
    read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  })
  # print(list_of_counts_data)
  return(list_of_counts_data)
}
read_common_genes <- function(file_path) {
  common_genes <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  return(common_genes$V1)
}
merge_counts_data <- function(list_of_counts_data, common_genes) {
  merged_data <- data.frame(GeneId = common_genes)
  for(counts_data in list_of_counts_data) {
    filtered_data <- subset(counts_data, GeneId %in% common_genes)
    merged_data <- merge(merged_data, filtered_data, by = "GeneId")
  }
  return(merged_data)
}
directory <- "Results"
pattern <- "*_Keep_Genes_Expression\\.txt$"
common_genes_file <- "Results/AT_common_significant_genes.txt"
list_of_counts_data <- read_counts_data(directory, pattern)
common_genes <- read_common_genes(common_genes_file)
merged_counts_data <- merge_counts_data(list_of_counts_data, common_genes)
merged_counts_data_name<- merged_counts_data[,1]
merged_counts_data <- merged_counts_data[,-1]
merged_counts_data <- log1p(merged_counts_data)		#ln(x+1)变换，尽力消除技术等差距，方便后期使用。
merged_counts_data$GeneId<- merged_counts_data_name
merged_counts_data <- merged_counts_data[, c("GeneId", setdiff(names(merged_counts_data), "GeneId"))]
output_file_path <- "Results/AT_merged_counts_data.txt"
write.table(merged_counts_data, file = output_file_path, sep = "\t", row.names = FALSE, quote = FALSE)

