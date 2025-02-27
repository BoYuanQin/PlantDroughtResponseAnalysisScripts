# R/wgcna_analysis.R
library(WGCNA)
library(flashClust)

# 修改后的类定义（关键修改点）
WeightedGeneCoexpressionNetworkAnalysis <- setRefClass(
  "WeightedGeneCoexpressionNetworkAnalysis",
  methods = list(
    plotSoftThreshold = function(sft) {
      output_file <- file.path(results_paths$WGCNA, "SoftThresholdPlot.pdf")
      create_dir_if_needed(dirname(output_file))
      pdf(output_file, width = 8, height = 6)
      # ...原有绘图代码...
    },
    
    buildNetwork = function(power) {
      # ...原有代码...
      
      # 修改所有输出路径
      write.table(TOM, file.path(results_paths$WGCNA, "TOM_matrix.txt"), 
                 sep = "\t", row.names = TRUE, col.names = TRUE)
      # ...其余输出路径修改...
    }
  )
)
