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
    # 加载数据
    loadData = function(filename) {
      if (!file.exists(filename)) stop("File not found!")
      exprMatrix <<- as.matrix(read.table(filename, header = TRUE, row.names = 1, sep = "\t"))
      exprMatrix <<- t(exprMatrix)
	  rn <- rownames(exprMatrix)
	  exprMatrix <<- apply(exprMatrix, 2, as.numeric)
	  rownames(exprMatrix) <<- rn
	  print(head(exprMatrix,n=10))
    }, 
	# 数据预处理，去除NA值，获取样本数量和设置网络类型
    preprocessData = function() {
      exprMatrix <<- na.omit(exprMatrix)
      exprMatrix <<- unique(exprMatrix)
      nSamples <<- nrow(exprMatrix)
      type <<- "unsigned"
    },
	# 选择软阈值，用于建立权重网络
    chooseSoftPower = function() {
      powers = c(c(1:10), seq(from = 12, to=30, by=2))  
      sft <- pickSoftThreshold(exprMatrix, powerVector = powers, networkType=type, verbose = 5)
	  plotSoftThreshold(sft)
	  sft_data <- data.frame(
		Power = powers,
		sft$fitIndices,
		PowerEstimate = rep(sft$powerEstimate, length(powers))
	  )
	  write.table(sft_data, file.path(results_paths$WGCNA, "sft_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      return(sft$powerEstimate) 
    },
	# 绘制软阈值图
	plotSoftThreshold = function(sft) {
      output_file <- file.path(results_paths$WGCNA, "SoftThresholdPlot.pdf")
      create_dir_if_needed(dirname(output_file))
      pdf(output_file, width = 8, height = 6)
      # pdf("PPI_WGCNA_Intersection_Genes/WGCNA/SoftThresholdPlot.pdf", width = 8, height = 6)
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
	# 建立基因表达网络，计算邻接矩阵和拓扑重叠矩阵(TOM)
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
      # write.table(TOM, file="PPI_WGCNA_Intersection_Genes/WGCNA/TOM_matrix.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
      write.table(TOM, file.path(results_paths$WGCNA, "TOM_matrix.txt"), 
           sep = "\t", row.names = TRUE, col.names = TRUE)
      dissTOM <- 1-TOM
      geneTree <<- flashClust(as.dist(dissTOM), method="average")
	  #保存geneTree
	  if (requireNamespace("ape", quietly = TRUE)) {
		library(ape)
		ape::write.tree(as.phylo(geneTree),file.path(results_paths$WGCNA, "geneTree.txt") )
	  } else {
		dendro <- as.dendrogram(geneTree)
		write.table(dendro, file.path(results_paths$WGCNA, "geneTree.txt"), sep = "\t", quote = FALSE)
	  }
      dynamicMods <<- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 10)	#就是每个基因的类别
	  # print(head(dynamicMods,10))
      moduleLabels <<- labels2colors(dynamicMods)
	  #保存moduleLables
	  write.table(moduleLabels,  file.path(results_paths$WGCNA, "moduleLabels.txt"), sep = "\t", col.names = NA, quote = FALSE)
      return(list(adjacencyMatrix, TOM, moduleLabels))
    },
	# 可视化构建的基因网络
    visualizeNetwork = function() {
      pdf(file.path(results_paths$WGCNA, "Original_Cluster_Tree.pdf"), width = 8/2.5, height = 6/2.5)
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
      write.table(new_data, file.path(results_paths$WGCNA, "Trait_Data.txt"), sep = '\t', row.names = FALSE, quote = FALSE)
    },
	# 加载和保存表型数据
    loadTraitData = function(filename) {
      if (!file.exists(filename)) stop("Trait data file not found!")
		traitData <<- read.table(filename, header = TRUE, row.names = 1, sep = "\t")
      if (!identical(rownames(traitData), rownames(exprMatrix))) {
        traitData <<- traitData[rownames(exprMatrix), , drop = FALSE]
      }
    },
	# 分析模块与表型特征之间的关系
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
	  write.table(.self$correlationAndPvalue, file.path(results_paths$WGCNA, "Module_Trait_Relationships.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    },
	# 选择关键基因
    selectKeyGenes = function() {
      if (!is.data.frame(.self$correlationAndPvalue)) {
		.self$correlationAndPvalue <- as.data.frame(.self$correlationAndPvalue)
	  }
	  print(correlationAndPvalue)
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
	# 保存结果
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


