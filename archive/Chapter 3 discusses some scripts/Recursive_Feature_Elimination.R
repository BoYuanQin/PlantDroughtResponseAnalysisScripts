library(caret)
library(ggplot2)
Recursive_Feature_Elimination <- setRefClass(
  "Recursive_Feature_Elimination",
  fields = list(
    X = "matrix",
    Y = "factor"
  ),
  methods = list(
    runRFE = function(Funcs) {
      number <- 5
      subsetSizes <- c(2:ncol(.self$X))
      set.seed(123)
      seeds <- vector(mode = "list", length = number + 1)
	  for(i in 1:number) {
	  	seeds[[i]] <- sample.int(1000, length(subsetSizes))
	  }
	  seeds[[number + 1]] <- sample.int(1000, 1)
      control <- rfeControl(functions = Funcs, 
							method = "cv", 
							number = number, 
							returnResamp = "all",
							allowParallel=TRUE,
                      	    seeds = seeds)
      results <- rfe( .self$X, .self$Y, sizes = subsetSizes, 
					  rfeControl = control, metric ="Accuracy",
					  maximize = TRUE)
      return(results)
    },
    plotAccuracy = function(rfeResults) {
      performanceData <- rfeResults$results
      max_Accuracy_index <- which.max(performanceData$Accuracy)
      max_Accuracy <- performanceData[max_Accuracy_index, ]
      ggplot(performanceData, aes(x = Variables, y = Accuracy)) + geom_line() + 
				geom_point(data = max_Accuracy, aes(x = Variables, y = Accuracy), color = "red", size = 4)+ 
				xlab("Number of Features") + ylab("Accuracy")     
    }
  )
)
svmFuncs <- caretFuncs
svmFuncs$fit <-function(x, y, first, last, ...) {
	train(as.matrix(x), y, method = "svmRadial", ...)	#svmRadial svmLinear
}
rfFuncs <- caretFuncs
rfFuncs$fit <-function(x, y, first, last, ...) {
    train(as.matrix(x), y, method = "rf", ...)
}
glmnetFuncs <- caretFuncs
glmnetFuncs$fit <-function(x, y, first, last, ...) {
    train(as.matrix(x), y, method = "glmnet", ...)
}
xgbTreeFuncs <- caretFuncs
xgbTreeFuncs$fit <-function(x, y, first, last, ...) {
    train(as.matrix(x), y, method = "xgbTree", ...)
}

GetIntersectionGenesData <- function(PPI_Key_Path, WGCNA_Key_Path, GeneExpressionPath) {
	PPI_Key_Genes <- read.table(PPI_Key_Path, header = FALSE,sep="\t",stringsAsFactors=FALSE, quote = "\"")[, 1]
	WGCNA_Key_Genes <- read.table(WGCNA_Key_Path, header = FALSE,sep="\t",stringsAsFactors=FALSE, quote = "\"")[, 1]
	keyGenes <- intersect(PPI_Key_Genes, WGCNA_Key_Genes)
	GeneExpression <- read.table(GeneExpressionPath, header = TRUE,sep="\t",stringsAsFactors=FALSE, quote = "\"")
	colnames(GeneExpression)[1] <- "GeneID"
	keyGenesExpression <- t(GeneExpression[GeneExpression$GeneID %in% keyGenes, ])
	colnames(keyGenesExpression) <- keyGenesExpression[1, ]
	keyGenesExpression <- keyGenesExpression[-1, ]
	# print(keyGenesExpression)
	write.table(keyGenesExpression, file = "Results/KeyGenesExpression.txt", sep = "\t", quote = FALSE, row.names = TRUE)	
	Y <- ifelse(grepl("control", rownames(keyGenesExpression)), "A", "B")
	rownames(keyGenesExpression) <- NULL
	X <- keyGenesExpression
	return(list(X = X, Y = Y))
}
TrainData <- GetIntersectionGenesData("Data/RFE-AT/ppi_key_genes.txt", "Data/RFE-AT/wgcna_KeyGenes.txt", "Data/RFE-AT/merged_counts_data.txt")
# X <- TrainData$X[, 1:10]
X <- TrainData$X
X <- apply(X, 2, as.numeric)
# print(X)
Y <- TrainData$Y
Y <- factor(Y)
# print(Y)


runAndSaveRFE <- function(RfeInstance, funcs, prefix) {
	#开始运行
	RfeInstance <- Recursive_Feature_Elimination$new()
	RfeInstance$X <- X
	RfeInstance$Y <- Y
	# warnings()
    # 运行 RFE
    rfe_results <- RfeInstance$runRFE(funcs)    
    # 保存 RFE 结果汇总
    rfe_summary <- rfe_results$results
    summary_filename <- paste(prefix, "summary.txt", sep = "_")
    write.table(rfe_summary, summary_filename, row.names = FALSE, sep = "\t", quote = FALSE)	    
    # 保存最佳性能对应的特征
    optimal_features <- predictors(rfe_results, rfe_results$optsize)
    features_filename <- paste(prefix, "optimal_features.txt", sep = "_")
	write.table(as.data.frame(optimal_features), features_filename, row.names = FALSE, sep = "\t", quote = FALSE)
    # 绘制准确率的图并保存
    accuracy_plot <- RfeInstance$plotAccuracy(rfe_results)
    plot_filename <- paste(prefix, "Accuracy_Plot.pdf", sep = "_")
    # ggsave(plot_filename, accuracy_plot, width = 8, height = 6)
	ggsave(plot_filename, plot = accuracy_plot, device = "pdf", width = 8, height = 6)
	return(optimal_features)
}
#此并行和函数内部并行冲突，不能执行并行。
# library(parallel)
# tasks <- list(
  # rf = function() runAndSaveRFE(RfeInstance, rfFuncs, "Results/rfFuncs"),
  # svm = function() runAndSaveRFE(RfeInstance, svmFuncs, "Results/svmFuncs"),
  # glmnet = function() runAndSaveRFE(RfeInstance, glmnetFuncs, "Results/glmnetFuncs")
# )
# results <- mclapply(tasks, mc.cores = detectCores()-1)
# rf_optimal_features <- results$rf
# svm_optimal_features <- results$svm
# glmnet_optimal_features <- results$glmnet

rf_optimal_features <- runAndSaveRFE(RfeInstance, rfFuncs, "Results/rfFuncs")
svm_optimal_features <- runAndSaveRFE(RfeInstance, svmFuncs, "Results/svmFuncs")
glmnet_optimal_features <- runAndSaveRFE(RfeInstance, glmnetFuncs, "Results/glmnetFuncs")
common_optimal_features <- Reduce(intersect, list(rf_optimal_features, svm_optimal_features, glmnet_optimal_features))
write.table(common_optimal_features ,file = "Results/common_optimal_features.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# runAndSaveRFE(RfeInstance, xgbTreeFuncs, "Results/xgbTreeFuncs")
