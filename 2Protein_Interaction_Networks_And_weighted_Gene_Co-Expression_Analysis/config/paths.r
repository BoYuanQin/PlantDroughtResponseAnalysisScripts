# config/paths.R
input_paths <- list(
  gene_list = "input/gene_list/common_significant_genes.txt",
  expression_data = "input/expression_data/deg_data_counts.tsv",
  string_db = "input/STRINGdb"  # 需预先下载STRING数据库
)

results_paths <- list(
  PPI = "results/PPI",
  WGCNA = "results/WGCNA",
  INTERSECTION = "results/Intersection"
)
