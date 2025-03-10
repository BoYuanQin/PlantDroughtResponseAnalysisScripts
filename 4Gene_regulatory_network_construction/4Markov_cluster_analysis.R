# conda activate genesinference
GeneNetworkAnalysis <- setRefClass(
	"GeneNetworkAnalysis",
	fields = list(
		Network_GENIE3_significant = "data.frame"
	),
	methods = list(
		initialize = function() {},
		load_network_from_file = function(file_path) {
		# 从文件中加载基因网络，参数 file_path 指定文件路径。
			.self$Network_GENIE3_significant <- read.table(file_path, header = TRUE, sep = "\t")
		},
		# 执行网络聚类和基因本体(GO)富集分析，输出结果到指定目录。
        # 参数包括输出目录、生物体(默认为玉米)、最小聚类大小、MCL算法的参数等。
		perform_network_clustering_and_go_analysis = function(output_dir, organism = "zmays", min_cluster_size=10, mcl_I = 2.5, mcl_te = 20, min_go_term_size = 10, max_p_value = 0.05) {
			library(tidyverse)
			library(stringr)
			library(gprofiler2)

			mcl_path <- "mcl"
			mcl_output_file <- paste0(output_dir, "/mcl_cluster_I", mcl_I, ".txt")
			system(paste(mcl_path, "OutputData/AT/GENIE3NetworkAll/Network_GENIE3_significant.txt", "--abc", paste("-I", mcl_I, "-te", mcl_te, "-o", mcl_output_file)))
			clust <- read_lines(file = mcl_output_file) %>% str_split(., "\t")
			clust <- clust[lengths(clust) > min_cluster_size]
			cleaned_genes <- lapply(clust, function(x) gsub("\"", "", x))
			for (i in seq_along(cleaned_genes)) {
				go_enrich <- gost(query = unlist(cleaned_genes[[i]]), organism = organism, user_threshold = max_p_value, correction_method = "g_SCS", domain_scope = "annotated")
				if (is.null(go_enrich$result) || nrow(go_enrich$result) == 0) {
					next
				}
				go_enrich_result <- go_enrich$result %>% arrange(p_value)
				go_enrich$result <- go_enrich_result

				plot_file <- paste0(output_dir, "/module_", i, "_go_enrich_plot.pdf")
				pdf(plot_file, width = 10, height = 6)
				p <- gostplot(go_enrich, capped = TRUE, interactive = FALSE)
				pp <- publish_gostplot(p, highlight_terms = go_enrich$result[c(1:5),], 
										width = NA, height = NA, filename = NULL )
				print(pp)
				dev.off()
				if (length(go_enrich_result) > 0) {
					go_output_file <- paste0(output_dir, "/module_", i, "_", "go_enrich.txt")
					write_tsv(go_enrich_result, file = go_output_file)                    
				}
			}
		},
		# 计算转录因子的度中心性，并绘制条形图。
        # 参数 output_file_barplot 指定输出文件路径。
		plot_tf_degree_centrality = function(output_file_barplot) {
			library(tidyverse)

			cal_deg = function(link_list) {
				link_list %>% select(-importance) %>% group_by(TF) %>% summarise(n = n()) %>% arrange(desc(n))
			}
		
			tf_degree_poisson <- cal_deg(.self$Network_GENIE3_significant)
			write.table(tf_degree_poisson, file = "OutputData/AT/EnrichmentResults/TFs_Degree.txt", sep = "\t", quote = FALSE, row.names = TRUE)
			lambda <- mean(tf_degree_poisson$n)
			tf_degree_poisson$poisson_prob <- dpois(tf_degree_poisson$n, lambda)
			write.table(tf_degree_poisson, file = "OutputData/AT/EnrichmentResults/tf_degree_poisson.txt", sep = "\t", quote = FALSE, row.names = TRUE)
			threshold <- 0.05
			key_tfs <- tf_degree_poisson[tf_degree_poisson$poisson_prob < threshold, ]
			write.table(key_tfs, file = "OutputData/AT/EnrichmentResults/key_tfs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
			pdf(output_file_barplot, width = 10, height = 6)
			p_tissue <- ggplot(data = tf_degree_poisson) +
						geom_col(mapping = aes(x = reorder(TF, -n), y = n)) +
						theme(plot.title = element_text(hjust = 0.5)) +
						labs(x = "Transcription Factors", y = "Number of targets", title = "TF Degree Centrality")
			print(p_tissue)
			dev.off()
		},
		# 读取并合并文件夹中的文件，去除重复项，并输出到统一文件。
        # 参数 directory 指定文件夹路径。
		read_and_merge_files = function(directory) {
			library(tidyverse)
			file_list <- list.files(directory, pattern = "^module_\\d+_go_enrich\\.txt$", full.names = TRUE)
			merged_data <- tibble()  
			for (file_name in file_list) {
				# total_lines <- count.fields(file_name, skip = 1)
				# lines_to_read <- max(1, ceiling(total_lines * 0.1))
				# # data <- read_tsv(file_name, n_max = lines_to_read)  
				data <- read_tsv(file_name)  
				data <- data %>% mutate(FileName = basename(file_name))
				merged_data <- bind_rows(merged_data, data)  
			}
			merged_data_unique <- merged_data %>% distinct(term_name, .keep_all = TRUE)
			write_tsv(merged_data_unique, "OutputData/AT/EnrichmentResults/merged_modules_enrichment_results_all.txt")
		}
	)
)
set.seed(123)
net <- GeneNetworkAnalysis$new()
net$load_network_from_file("OutputData/AT/GENIE3NetworkAll/Network_GENIE3_significant.txt")
net$perform_network_clustering_and_go_analysis(output_dir="OutputData/AT/EnrichmentResults", organism = "athaliana") #athaliana zmays osativa
net$plot_tf_degree_centrality(output_file_barplot="OutputData/AT/EnrichmentResults/barplot.pdf")
net$read_and_merge_files("OutputData/AT/EnrichmentResults")

