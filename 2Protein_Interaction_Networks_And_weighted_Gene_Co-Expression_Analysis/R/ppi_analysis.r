# conda acitvate newproteinanalysis
# 加载必要的库
library(STRINGdb)
library(igraph)
library(MCL)
library(ProNet)
library(parallel)
options(warn=-1)

create_dir_if_needed <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# 将 save_subnetworks 函数移出类定义，作为独立的函数

save_subnetworks <- function(main_graph, key_modules_names, format = "graphml", path) {
  extract_subnetwork <- function(main_graph, vertices) {
    induced_subgraph(main_graph, vids = V(main_graph)[name %in% vertices])
  }
  
  create_dir_if_needed(path)  
  
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


# 定义 WgcnaNetworkAnalyzer 类
WgcnaNetworkAnalyzer <- setRefClass(
  "WgcnaNetworkAnalyzer",

  fields = list(
    geneIDs = "ANY",          
    proteinIDs = "character",
    string_db = "ANY",               
    PPI_subnetwork = "ANY",         
    network_modules = "list",       
    key_modules = "list",           
    key_genes = "list",        
    dat_map = "ANY"
  ),

  methods = list(
    read_geneIDs = function(file_path) {
      if(!file.exists(file_path)) {
        stop("File does not exist!")
      }
      gene_data <- read.table(file_path, header=TRUE, stringsAsFactors = FALSE)
      .self$geneIDs <- gene_data[,1]
      if(length(geneIDs) == 0) {
        stop("No gene IDs found!")
      }
      cat("Gene IDs read successfully.\n")
    },

    get_string_data = function(local_path) {
      if(file.exists(local_path)) {
        string_db <<- STRINGdb$new(version="12.0", species=3702, score_threshold=900, network_type="full", input_directory=local_path)
        cat("STRING data loaded from local directory successfully.\n")
      } else {
        string_db <<- STRINGdb$new(version="12.0", species=3702, score_threshold=900, network_type="full", input_directory=local_path)
        cat("STRING data saved to local directory.\n")
      }
    },

    map_gene_to_protein = function() {
      dat_map <<- string_db$map(my_data_frame=data.frame(gene_id=geneIDs), 
                                my_data_frame_id_col_names="gene_id", 
                                removeUnmappedRows = TRUE)
      if(nrow(dat_map) == 0) {
        stop("No STRING IDs mapped!")
      }
      proteinIDs <<- dat_map$STRING_id
      cat("Gene IDs mapped to STRING IDs successfully using STRINGdb.\n")
    },

    construct_subnetwork = function() {
      if(length(proteinIDs) == 0) {
        stop("No protein IDs provided!")
      }
      PPI_network <- string_db$get_subnetwork(proteinIDs)
      if(gorder(PPI_network) == 0) {
        stop("No PPI edges found for the provided protein IDs!")
      }
      PPI_subnetwork <<- PPI_network
      PPI_subnetwork <<- simplify(PPI_subnetwork, remove.multiple = TRUE, remove.loops = TRUE)
      cat("PPI subnetwork constructed successfully using STRINGdb functions.\n")
    },

    cluster_network_iGraph = function() {
      if (is.null(PPI_subnetwork) || gorder(PPI_subnetwork) == 0) {
        stop("PPI subnetwork is not available or empty!")
      }

      clustering_functions <- list(
        fastgreedy = function(graph) {
          cluster_result <- cluster_fast_greedy(graph)
          membership_vec <- membership(cluster_result)
          cluster_list <- split(V(graph)$name, membership_vec)
          output <- capture.output({
            cat("fastgreedy clusters (first 5):\n")
            print(cluster_list[1:min(5, length(cluster_list))])
          })
          return(list(clusters = cluster_list, output = output))
        },
        MCL = function(graph) {
          adj_matrix <- as.matrix(as_adjacency_matrix(graph))
		      # write.csv(adj_matrix, file.path(results_paths$PPI, "adj_matrix.csv"))
          mcl_result <- MCL::mcl(adj_matrix, addLoops = TRUE, inflation = 2.5, max.iter =100)
		      # write.csv(mcl_result, "mcl_result.csv")
          # write.csv(mcl_result, file.path(results_paths$PPI, "mcl_result.csv"))
          output <- capture.output({
            cat("MCL mcl_result (first 5):\n")
            print(mcl_result[1:min(5, length(mcl_result))])
          })
          cluster_assignments <- mcl_result$Cluster
          names(cluster_assignments) <- rownames(adj_matrix)
          cluster_list <- split(names(cluster_assignments), cluster_assignments)
          return(list(clusters = cluster_list, output = output))
        },
        MCODE = function(graph) {
          mcode_result <- mcode(graph, vwp=0.05, haircut=TRUE, fluff=FALSE, fdt=0.8, loops=FALSE)
          output <- capture.output({
            cat("MCODE mcode_result (first 5):\n")
            print(mcode_result[1:min(5, length(mcode_result))])
          })
          cluster_indices <- mcode_result$COMPLEX
          if (length(cluster_indices) == 0) {
            cluster_list <- list()
          } else {
            cluster_list <- lapply(cluster_indices, function(module) {
              V(graph)$name[module]
            })
          }
          return(list(clusters = cluster_list, output = output))
        }
      )

      algorithms <- names(clustering_functions)
      num_cores <- min(length(algorithms), detectCores())
      cat("Using", num_cores, "cores for parallel execution.\n")

      cl <- makeCluster(num_cores)
      clusterExport(cl, varlist = c("PPI_subnetwork", "clustering_functions", "save_subnetworks", "create_dir_if_needed",  "results_paths"), envir = environment())
      clusterEvalQ(cl, {
        library(igraph)
        library(MCL)
        library(ProNet)
      })

      results_list <- parLapply(cl, algorithms, function(algo) {
        cat("Running", algo, "clustering algorithm.\n")
        result <- clustering_functions[[algo]](PPI_subnetwork)
        cluster_list <- result$clusters
        output <- result$output

        # 移除节点数小于等于10的模块
        cluster_list <- cluster_list[sapply(cluster_list, length) > 10]

        # 检查 clusters 是否为空
        if (length(cluster_list) == 0) {
          output <- c(output, paste("No clusters found for algorithm:", algo))
          return(list(algo = algo, clusters = cluster_list, output = output))
        } else {
          output <- c(output, paste("Number of clusters found by", algo, ":", length(cluster_list)))
        }

        # 保存子网络和模块大小
        # dir.create(paste0("./", algo), showWarnings = FALSE)
        dir_path <- file.path(results_paths$PPI, algo)
        dir.create(dir_path, showWarnings = FALSE, recursive = TRUE) 
        save_subnetworks(main_graph=PPI_subnetwork, key_modules_names=cluster_list, format = "graphml", path=file.path(results_paths$PPI, algo))
        module_sizes <- sapply(cluster_list, length)
        modules_df <- data.frame(
          module_name = paste0("module_", seq_along(module_sizes)),
          size = module_sizes
        )
        write.table(modules_df, file =file.path(results_paths$PPI,  paste0(algo, "_module_sizes.txt")), sep = "\t", row.names = FALSE)

        return(list(algo = algo, clusters = cluster_list, output = output))
      })

      stopCluster(cl)

      # 在主进程中打印子进程的输出
      for (res in results_list) {
        cat(paste("Algorithm:", res$algo, "\n"))
        cat(paste(res$output, collapse = "\n"), "\n")
      }

      # 组织结果
      key_modules_temp <- list()
      for (res in results_list) {
        key_modules_temp[[res$algo]] <- list(clusters = res$clusters)
      }
      key_modules <<- key_modules_temp
    },

    get_key_genes = function() {
      if(is.null(key_modules) || length(key_modules) == 0) {
        stop("No key modules found!")
      }
      key_genes_list <- list()
      for (algo in names(key_modules)) {
        cluster_list <- key_modules[[algo]]$clusters
        if (length(cluster_list) == 0) {
          cat("No clusters to process for algorithm:", algo, "\n")
          next
        }
        genes_with_module <- data.frame(gene = character(0), module = character(0))
        for (i in seq_along(cluster_list)) {
          module_genes <- cluster_list[[i]]
          mapped_genes <- dat_map$gene_id[dat_map$STRING_id %in% module_genes]
          genes_with_module <- rbind(genes_with_module, data.frame(gene = na.omit(mapped_genes), module = paste0(algo, "_module_", i)))
        }
        geneIDs_upper <- toupper(geneIDs)
        genes_with_module$gene <- toupper(genes_with_module$gene)
        genes_with_module <- genes_with_module[genes_with_module$gene %in% geneIDs_upper, ]
        key_genes_list[[algo]] <- genes_with_module
      }
      key_genes <<- key_genes_list
      cat("Key genes obtained successfully using dat_map.\n")
    },

    save_key_genes = function(file_path) {
      create_dir_if_needed(dirname(file_path))  
      
      lapply(names(key_genes), function(algo) {
        output_file <- file.path(
          dirname(file_path), 
          paste0(algo, "_", basename(file_path))
        )
        write.table(key_genes[[algo]], output_file, 
                   sep = "\t", row.names = FALSE, col.names = TRUE)
      })
      cat("Key genes saved to:", dirname(file_path), "\n")
    }
  )
)

