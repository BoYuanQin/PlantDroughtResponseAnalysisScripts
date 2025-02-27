# R/ppi_analysis.R
library(STRINGdb)
library(igraph)
library(MCL)
library(ProNet)
library(parallel)

# 新增：动态路径处理函数
create_dir_if_needed <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# 修改后的保存子网络函数
save_subnetworks <- function(main_graph, key_modules_names, format = "graphml", path) {
  extract_subnetwork <- function(main_graph, vertices) {
    induced_subgraph(main_graph, vids = V(main_graph)[name %in% vertices])
  }
  
  create_dir_if_needed(path)  # 新增目录创建检查
  
  lapply(seq_along(key_modules_names), function(i) {
    module <- key_modules_names[[i]]
    subgraph <- extract_subnetwork(main_graph, module)
    filename <- file.path(path, paste0("module_", i, ".", format))
    
    if (format == "graphml") {
      V(subgraph)$label <- V(subgraph)$name
      write_graph(subgraph, filename, "graphml")
    } else if (format == "sif") {
      write_graph(subgraph, filename, "edgelist")
    }
  })
}

# 修改后的类定义（关键修改点）
WgcnaNetworkAnalyzer <- setRefClass("WgcnaNetworkAnalyzer",
  methods = list(
    save_key_genes = function(file_path) {
      create_dir_if_needed(dirname(file_path))  # 新增目录检查
      
      lapply(names(key_genes), function(algo) {
        output_file <- file.path(
          dirname(file_path), 
          paste0(algo, "_", basename(file_path))
        )
        write.table(key_genes[[algo]], output_file, 
                   sep = "\t", row.names = FALSE, col.names = TRUE)
      })
      cat("Key genes saved to:", dirname(file_path), "\n")
    },
    
    cluster_network_iGraph = function() {
      # ...原有代码...
      
      # 修改输出路径
      lapply(results_list, function(res) {
        algo_dir <- file.path(results_paths$PPI, res$algo)
        create_dir_if_needed(algo_dir)
        
        save_subnetworks(PPI_subnetwork, res$clusters, 
                        "graphml", algo_dir)
        # ...其余保存逻辑...
      })
    }
  )
)
