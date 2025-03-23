# Load necessary libraries
library(igraph)
library(poweRlaw)

# Function to read data and build the graph
build_graph <- function(file_path) {
  # Read the data from the specified file path
  data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Build the graph using the 'TF' and 'target' columns
  g <- graph_from_data_frame(d = data[, c("TF", "target")], directed = TRUE)
  return(g)
}

# Function to assess scale-free property
assess_scale_free <- function(g, output_prefix) {
  # Calculate degrees of all nodes
  degrees <- degree(g, mode = "all")
  
  # Remove nodes with zero degree for power-law fitting
  degrees_nonzero <- degrees[degrees > 0]
  
  # Fit power-law model using poweRlaw package
  m_pl <- displ$new(degrees_nonzero)
  est <- estimate_xmin(m_pl)
  m_pl$setXmin(est)
  
  # Plot the degree distribution with power-law fit and save the plot
  # png(filename = paste0(output_prefix, "_degree_distribution.png"))
  pdf(file = paste0(output_prefix, "_degree_distribution.pdf"))
  plot(m_pl, main = "Degree Distribution with Power-law Fit", xlab = "Degree", ylab = "Frequency")
  lines(m_pl, col = "red")
  dev.off()
  
  # Goodness-of-fit test
  num_cores <- parallel::detectCores() - 1
  gof <- bootstrap_p(m_pl, threads = num_cores, no_of_sims = 1000)
  
  # Save results to a list
  results <- list(
    alpha = m_pl$pars,
    xmin = m_pl$xmin,
    p_value = gof$p
  )
  
  # Save results to a text file
  results_df <- data.frame(
    Alpha = round(results$alpha, 3),
    Xmin = results$xmin,
    P_value = round(results$p_value, 3)
  )
  write.table(results_df, file = paste0(output_prefix, "_scale_free_results.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  
  return(results)
}

# Function to assess small-world property
assess_small_world <- function(g, output_prefix) {
  # Calculate average clustering coefficient
  avg_clustering <- transitivity(g, type = "average")
  
  # Calculate average shortest path length
  if (is.connected(g)) {
    avg_path_length <- average.path.length(g, directed = TRUE)
  } else {
    avg_path_length <- mean_distance(g, directed = TRUE, unconnected = TRUE)
  }
  
  # Generate a random graph with the same degree sequence
  deg_seq_in <- degree(g, mode = "in")
  deg_seq_out <- degree(g, mode = "out")
  g_random <- degree.sequence.game(out.deg = deg_seq_out, in.deg = deg_seq_in, method = "simple")
  
  # Calculate clustering coefficient and average path length for the random graph
  avg_clustering_random <- transitivity(g_random, type = "average")
  if (is.connected(g_random)) {
    avg_path_length_random <- average.path.length(g_random, directed = TRUE)
  } else {
    avg_path_length_random <- mean_distance(g_random, directed = TRUE, unconnected = TRUE)
  }
  
  # Save results to a list
  results <- list(
    avg_clustering = avg_clustering,
    avg_path_length = avg_path_length,
    avg_clustering_random = avg_clustering_random,
    avg_path_length_random = avg_path_length_random
  )
  
  # Save results to a text file
  results_df <- data.frame(
    Metric = c("Average Clustering Coefficient", "Average Path Length"),
    Original_Network = round(c(avg_clustering, avg_path_length), 3),
    Random_Network = round(c(avg_clustering_random, avg_path_length_random), 3)
  )
  write.table(results_df, file = paste0(output_prefix, "_small_world_results.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  
  return(results)
}

# Main function to run the analysis
analyze_network <- function(file_path, output_prefix) {
  # Build the graph
  g <- build_graph(file_path)
  
  # Assess scale-free property
  scale_free_results <- assess_scale_free(g, output_prefix)
  
  # Assess small-world property
  small_world_results <- assess_small_world(g, output_prefix)
  
  # Save combined results to a list
  all_results <- list(
    scale_free = scale_free_results,
    small_world = small_world_results
  )
  
  # Optionally, save all results to an RData file
  save(all_results, file = paste0(output_prefix, "_all_results.RData"))
  
  return(all_results)
}


analyze_network("Network_GENIE3_significant.txt", "network_analysis_output")
