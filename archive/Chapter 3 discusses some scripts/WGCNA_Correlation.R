library(WGCNA)
library(flashClust)

enableWGCNAThreads(nThreads = 3)

WeightedGeneCoexpressionNetworkAnalysis <- setRefClass(
  "WeightedGeneCoexpressionNetworkAnalysis",

  fields = list(
    exprMatrix = "matrix",  
    softPower = "numeric",   
    adjacencyMatrix = "matrix", 
    TOM = "matrix", 
    moduleLabels = "character", 
    nSamples = "numeric",
    type = "character",
    traitData = "data.frame", 
    moduleTraitPvalue = "ANY", 
	correlationAndPvalue= "data.frame", 
    geneTree= "ANY",
    dynamicMods= "ANY"),
  
  methods = list(
    loadData = function(filename) {
      if (!file.exists(filename)) stop("File not found!")
      exprMatrix <<- as.matrix(read.table(filename, header = TRUE, row.names = 1, sep = "\t"))
      exprMatrix <<- t(exprMatrix)
	  rn <- rownames(exprMatrix)
	  exprMatrix <<- apply(exprMatrix, 2, as.numeric)
	  rownames(exprMatrix) <<- rn
	  print(head(exprMatrix,n=10))
    }, 
    preprocessData = function() {
      exprMatrix <<- na.omit(exprMatrix)
      exprMatrix <<- unique(exprMatrix)
      nSamples <<- nrow(exprMatrix)
      type <<- "unsigned"
    },
    chooseSoftPower = function() {
      powers = c(c(1:10), seq(from = 12, to=30, by=2))  
      sft <- pickSoftThreshold(exprMatrix, powerVector = powers, networkType=type, verbose = 5)
	  plotSoftThreshold(sft)
	  sft_data <- data.frame(
		Power = powers,
		sft$fitIndices,
		PowerEstimate = rep(sft$powerEstimate, length(powers))
	  )
	  write.table(sft_data, file = "Results/sft_results.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      return(sft$powerEstimate) 
    },
	plotSoftThreshold = function(sft) {
      pdf("Results/SoftThresholdPlot.pdf", width = 8, height = 6)
      par(mfrow = c(1,2))
      powers = c(c(1:10), seq(from = 12, to=30, by=2))
      cex1 = 0.9
      plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
           main = paste("Scale independence"))
      text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
	   
      abline(h=0.80,col="red")
	  
      plot(sft$fitIndices[,1], sft$fitIndices[,5],
           xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
           main = paste("Mean connectivity"))
      text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
      dev.off()
    },
    buildNetwork = function(power) {
      if (is.na(power)){
        power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                ifelse(type == "unsigned", 22, 44))       
                )
                )
      }       
      softPower <<- power      
      print(paste("软阈值是：", softPower))
      adjacencyMatrix <<- adjacency(exprMatrix, power = softPower)
      adjacencyMatrix[is.na(adjacencyMatrix)] <<- 0
      TOM <<- TOMsimilarityFromExpr(exprMatrix, power = softPower)
      TOM[is.na(TOM)] <<- 0
      rownames(TOM) <<- colnames(exprMatrix)
      colnames(TOM) <<- colnames(exprMatrix)
      write.table(TOM, file="Results/TOM_matrix.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
      dissTOM <- 1-TOM
      geneTree <<- flashClust(as.dist(dissTOM), method="average")
	  #保存geneTree
	  if (requireNamespace("ape", quietly = TRUE)) {
		library(ape)
		ape::write.tree(as.phylo(geneTree), file = "Results/geneTree.txt")
	  } else {
		dendro <- as.dendrogram(geneTree)
		write.table(dendro, file = "Results/geneTree.txt", sep = "\t", quote = FALSE)
	  }
      dynamicMods <<- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 10)	#就是每个基因的类别
	  # print(head(dynamicMods,10))
      moduleLabels <<- labels2colors(dynamicMods)
	  #保存moduleLables
	  write.table(moduleLabels, file = "Results/moduleLabels.txt", sep = "\t", col.names = NA, quote = FALSE)
      return(list(adjacencyMatrix, TOM, moduleLabels))
    },
    visualizeNetwork = function() {
	  pdf("Results/Original_Cluster_Tree.pdf", width = 8/2.5, height = 6/2.5)
      plotDendroAndColors(geneTree, moduleLabels,"Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,cex.colorLabels = 0.6,cex.dendroLabels = 0.7)
      dev.off()
    },
    CreateAndSaveTreatmentData = function(file_path) {
      data <- read.table(file_path, header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)
      new_data <- data.frame(SamplesName = colnames(data))
	  conditions <- ifelse(grepl("control", new_data$SamplesName), "1", "0")
	  new_data$Treatment <- conditions
      # control_number <- as.integer(grepl("control", new_data$SamplesName))
      # new_data$Treatment <- 1 - control_number
      write.table(new_data, 'Results/Trait_Data.txt', sep = '\t', row.names = FALSE, quote = FALSE)
    },
    loadTraitData = function(filename) {
      if (!file.exists(filename)) stop("Trait data file not found!")
		traitData <<- read.table(filename, header = TRUE, row.names = 1, sep = "\t")
      if (!identical(rownames(traitData), rownames(exprMatrix))) {
        traitData <<- traitData[rownames(exprMatrix), , drop = FALSE]
      }
    },
    analyzeModuleTraitRelationship = function() {
      MEsData <- moduleEigengenes(exprMatrix, colors = moduleLabels,softPower=softPower)
	  MEs <- MEsData$eigengenes
      MEs = orderMEs(MEs)
	  print(MEs)
      moduleTraitCor <- cor(MEs, traitData, use = "p")
	  .self$moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
	  print(moduleTraitPvalue)
	  .self$correlationAndPvalue <- data.frame(Correlation = as.vector(moduleTraitCor), PValue = as.vector(moduleTraitPvalue))
	  rownames(.self$correlationAndPvalue) <- colnames(MEs)
	  write.table(.self$correlationAndPvalue, file = "Results/Module_Trait_Relationships_Cor.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    },
	analyzeModuleTraitRelationship_glm = function() {
      MEsData <- moduleEigengenes(exprMatrix, colors = moduleLabels,softPower=softPower)
	  MEs <- MEsData$eigengenes
      MEs = orderMEs(MEs)
	  # print(MEs)
      logistic_results <- list()
	  # pdf("Results/Model_Diagnostic_Plots.pdf", width = 8.27, height = 11.69) 
      for (me in colnames(MEs)) {
		# print(MEs[, me])
        model <- glm(as.factor(traitData$Treatment) ~ MEs[, me], family = binomial())
		pdf(paste("Results/Model_Diagnostic_Plots", me, ".pdf", sep="_"), width = 8.27, height = 11.69)
		par(mfrow=c(2, 2))  
		# 绘制模型诊断图
		plot(model, col='lightblue')  # 点的颜色设置为红色，按需调整
		# 添加一个自定义的标题
		title("Model Diagnostic Plots", col.main="black", font.main=2)  # 标题颜色和字体样式
		dev.off()
        summary_model <- summary(model)
        logistic_results[[me]] <- summary_model
      }      
      
	  coefficients_matrix <- sapply(logistic_results, function(x) coef(x)[2, "Estimate"])   
      pvalues_matrix <- sapply(logistic_results, function(x) coef(x)[2, "Pr(>|z|)"])
      combined_results <- data.frame(
		Estimate = coefficients_matrix,
		Pr = pvalues_matrix
	  )
	  write.table(combined_results, file = "Results/Module_Trait_Relationships_glm.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    },
    selectKeyGenes = function() {
      if (!is.data.frame(.self$correlationAndPvalue)) {
		.self$correlationAndPvalue <- as.data.frame(.self$correlationAndPvalue)
	  }
	  # print(correlationAndPvalue)
	  significantModules <- rownames(.self$correlationAndPvalue[.self$correlationAndPvalue$PValue < 0.05,])
	  
      keyGenes <- list()
      for (module in significantModules) {
        module <- sub("^ME", "", module)
        matchedIndices <- which(moduleLabels == module)
        genesInModule <- colnames(exprMatrix)[matchedIndices]
        keyGenes[[module]] <- genesInModule
      }      
      return(keyGenes)
    },    
    saveResults = function(filename) {
      keyGenes <- selectKeyGenes()
      con <- file(filename, "w")
      for (module in names(keyGenes)) {
        cat(paste0("Module: ", module, "\n"), file = con)
        write.table(keyGenes[[module]], file = con, sep = "\t", col.names =TRUE, row.names = FALSE, quote = FALSE, append = TRUE)
        cat("\n", file = con)
      }
      close(con)
    }    
  )
)

wgcna <- WeightedGeneCoexpressionNetworkAnalysis$new()
wgcna$loadData("Results/AT_merged_counts_data.txt")
wgcna$preprocessData()
power <- wgcna$chooseSoftPower()
network <- wgcna$buildNetwork(power)
wgcna$visualizeNetwork()
wgcna$CreateAndSaveTreatmentData("Results/AT_merged_counts_data.txt")
wgcna$loadTraitData("Results/Trait_Data.txt")
wgcna$analyzeModuleTraitRelationship()
wgcna$analyzeModuleTraitRelationship_glm()
wgcna$saveResults("Results/wgcna_KeyGenes.txt")


