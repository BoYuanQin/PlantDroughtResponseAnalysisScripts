# main.R
source("config/paths.R")
source("R/ppi_analysis.R")
source("R/wgcna_analysis.R")

# 创建输出目录
dir.create(results_paths$PPI, recursive = TRUE, showWarnings = FALSE)
dir.create(results_paths$WGCNA, recursive = TRUE, showWarnings = FALSE)
dir.create(results_paths$INTERSECTION, recursive = TRUE, showWarnings = FALSE)

# PPI分析流程
run_ppi_analysis <- function() {
  analyzer <- WgcnaNetworkAnalyzer$new()
  analyzer$read_geneIDs(input_paths$gene_list)
  analyzer$get_string_data(input_paths$string_db)
  analyzer$map_gene_to_protein()
  analyzer$construct_subnetwork()
  analyzer$cluster_network_iGraph()
  analyzer$get_key_genes()
  analyzer$save_key_genes(file.path(results_paths$PPI, "ppi_key_genes.txt"))
}

# WGCNA分析流程
run_wgcna_analysis <- function() {
  wgcna <- WeightedGeneCoexpressionNetworkAnalysis$new()
  wgcna$loadData(input_paths$expression_data)
  wgcna$preprocessData()
  power <- wgcna$chooseSoftPower()
  wgcna$buildNetwork(power)
  wgcna$visualizeNetwork()
  wgcna$CreateAndSaveTreatmentData(input_paths$expression_data)
  wgcna$loadTraitData(file.path(results_paths$WGCNA, "Trait_Data.txt"))
  wgcna$analyzeModuleTraitRelationship()
  wgcna$saveResults(file.path(results_paths$WGCNA, "wgcna_KeyGenes.txt"))
}

# 执行分析流程
if (!interactive()) {
  run_ppi_analysis()
  run_wgcna_analysis()
  
  # 调用Python脚本进行交集分析
  system(paste(
    "python3 python/find_intersection.py",
    "--wgcna", file.path(results_paths$WGCNA, "wgcna_KeyGenes.txt"),
    "--ppi_dir", results_paths$PPI,
    "--output_dir", results_paths$INTERSECTION
  ))
  
  message("Analysis completed! Results saved in results/ directory")
}
