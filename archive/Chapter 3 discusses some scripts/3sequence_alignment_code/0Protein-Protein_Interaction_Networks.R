# Load necessary libraries
library(STRINGdb)
library(igraph)
options(warn=-1)

WgcnaNetworkAnalyzer <- setRefClass(
  "WgcnaNetworkAnalyzer",
  
  fields = list(
    geneIDs = "ANY",          
    proteinIDs = "character",
    string_db = "ANY",               
    PPI_subnetwork = "ANY",         
    network_modules = "list",       
    key_modules = "list",           
    key_genes = "data.frame",        
    dat_map="ANY"
  ),
  methods = list(
    # 读取基因ID
    read_geneIDs = function(file_path) {
      if(!file.exists(file_path)) {
        stop("File does not exist!")
      }
      .self$geneIDs <- read.table(file_path, header=TRUE, stringsAsFactors = FALSE)[,1]
      # print(head(.self$geneIDs,n=10))
      if(length(geneIDs) == 0) {
        stop("No gene IDs found!")
      }
      cat("Gene IDs read successfully.\n")
    },
	# 从本地路径加载STRING数据库数据
    get_string_data = function(local_path) {
      if(file.exists(local_path)) {
        string_db <<- STRINGdb$new(version="12.0", species=4577, score_threshold=900, network_type="full", input_directory=local_path)
        cat("STRING data loaded from local RDS file using input_directory successfully.\n")
      } else {
        string_db <<- STRINGdb$new(version="12.0", species=4577, score_threshold=900, network_type="full", input_directory=local_path)	#ZM：4577 OS：39947 AT：3702	（物种代码）
        cat("STRING data saved to RDS file.\n")
      }
    },
	# 将基因ID映射到STRING蛋白质ID
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
	# 构建基于STRING ID的蛋白质相互作用子网络
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
	# 使用igraph进行网络聚类
    cluster_network_iGraph = function() {
      # 检查PPI子网络是否为空或不存在
	  if (is.null(PPI_subnetwork) || gorder(PPI_subnetwork) == 0) {
	  	stop("PPI subnetwork is not available or empty!")
	  }
	  # 执行聚类和评估
	  perform_clustering_and_evaluation <- function(graph, algo) {
	 	 # 处理spinglass算法时，选择最大组件
		 if (algo == "spinglass") {
	 	 components <- decompose.graph(graph)
	 	 graph <- components[[which.max(sapply(components, vcount))]]
	 	 }
	 	 cluster_result <- clustering_functions[[algo]](graph)
	 	 evaluation <- modularity(graph, membership = cluster_result$membership)
	 	 return(list(cluster = cluster_result, evaluation = evaluation))
	  }
	  # 定义聚类方法
	  clustering_functions <- list(
		  fastgreedy = cluster_fast_greedy,
		  walktrap = cluster_walktrap,
		  spinglass = cluster_spinglass,
		  edge_betweenness = cluster_edge_betweenness
	  )
	  # 使用并行计算进行多种聚类方法的评估
	  library(parallel)
	  results <- mclapply(names(clustering_functions), function(algo) {
	 	 perform_clustering_and_evaluation(PPI_subnetwork, algo)
	  }, mc.cores = detectCores())
	   # 评估各种聚类方法的效果，并选择最佳的聚类结果
	  evaluations <- sapply(results, function(x) x$evaluation)
	  best_clustering <- results[[which.max(evaluations)]]$cluster
	  # 打印聚类结果和子网络信息
	  print(best_clustering)
	  print(PPI_subnetwork)
	  # print(class(membership(best_clustering)))
	  # 保存聚类成员信息和顶点名
	  write.table(membership(best_clustering), file = "PPI_WGCNA_Intersection_Genes/PPI/best_clustering_check.txt", sep = "\t", quote = FALSE)
	  write.table(V(PPI_subnetwork)$name, file = "PPI_WGCNA_Intersection_Genes/PPI/PPI_subnetwork_Vertices_check.txt", sep = "\t", quote = FALSE)
	  #处理spinglass的情况 # 将顶点按聚类结果进行分类
	  V(PPI_subnetwork)$cluster <- NA
	  vertex_names <- V(PPI_subnetwork)$name
	  clusters <- sapply(vertex_names, function(v) {
        if (v %in% names(membership(best_clustering))) {
          return(membership(best_clustering)[v])
        } else {
          return(NA)      
        }})
	  V(PPI_subnetwork)$cluster <- clusters
	  # 移除未分类的顶点
	  PPI_subnetwork <<- delete.vertices(PPI_subnetwork, which(is.na(V(PPI_subnetwork)$cluster)))
	  key_modules <<- split(V(PPI_subnetwork)$name, V(PPI_subnetwork)$cluster)  
	  #PPI_subnetwork可以理解为稀疏邻接矩阵，列名是网络顶点名字，也是蛋白质的名字。
	  # 按聚类结果分割顶点，只保留大于10个顶点的聚类
	  key_modules <<- key_modules[sapply(key_modules, length) > 10]
	  #存储关键模块
	  save_subnetworks(main_graph=PPI_subnetwork, key_modules_names=key_modules, format = "graphml", path="PPI_WGCNA_Intersection_Genes/PPI")
	  #存储模块的大小
	  module_sizes <- sapply(key_modules, length)
	  modules_df <- data.frame(
	      module_name = names(module_sizes),
	      size = module_sizes
	  )
	  write.table(modules_df, file = "PPI_WGCNA_Intersection_Genes/PPI/module_sizes.txt", sep = "\t", row.names = FALSE)
    },    
	# 获取关键基因
    get_key_genes = function() {
      if(is.null(key_modules) || length(key_modules) == 0) {
        stop("No key modules found!")
      }
      genes_with_module <- data.frame(gene = character(0), module = integer(0))
      for (i in 1:length(key_modules)) {
        module_genes <- key_modules[[i]]
        mapped_genes <- dat_map$gene_id[dat_map$STRING_id %in% module_genes]
        genes_with_module <- rbind(genes_with_module, data.frame(gene = na.omit(mapped_genes), module = i))
      }
      geneIDs_upper <- toupper(geneIDs)
	  genes_with_module$gene <- toupper(genes_with_module$gene)
      key_genes <<- genes_with_module[genes_with_module$gene %in% geneIDs_upper, ]
      cat("Key genes obtained successfully using dat_map.\n")
    },  
	# 保存关键基因
    save_key_genes = function(file_path) {
      if(length(key_genes) == 0) {
        stop("No key genes found!")
      }
      write.table(key_genes, file_path, sep="\t",row.names = FALSE, col.names = FALSE)
      cat("Key genes saved successfully.\n")
    },
	# 保存子网络
    save_subnetworks = function(main_graph, key_modules_names, format = "graphml", path) {
		extract_subnetwork = function(main_graph, vertices) {
			subgraph <- induced_subgraph(main_graph, vids = V(main_graph)[name %in% vertices])
			return(subgraph)
		}
      for (i in 1:length(key_modules_names)) {
        module <- key_modules_names[[i]]
        subgraph <- extract_subnetwork(main_graph, module)
        filename <- paste0(path, "/module_", i, ".", format)
        if (format == "graphml") {
          write.graph(subgraph, file = filename, format = "graphml")
        } else if (format == "sif") {
          write.graph(subgraph, file = filename, format = "edgelist")
        }
      }
    }
  )
)

deg_analysis <- WgcnaNetworkAnalyzer$new() 
deg_analysis$read_geneIDs("Differentially_Expressed_Genes/common_significant_genes.txt")
deg_analysis$get_string_data("PPI_WGCNA_Intersection_Genes/PPI")
deg_analysis$map_gene_to_protein()    
deg_analysis$construct_subnetwork()    
deg_analysis$cluster_network_iGraph()
deg_analysis$get_key_genes()
deg_analysis$save_key_genes("PPI_WGCNA_Intersection_Genes/PPI/ppi_key_genes.txt")



