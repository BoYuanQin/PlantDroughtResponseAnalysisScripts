file_path <- "Results/玉米/raw_combined_data.txt"  
data <- read.table(file_path, header = TRUE)
col_names <- colnames(data)
labels <- ifelse(grepl("control", col_names, ignore.case = TRUE), "0", "1")
output <- data.frame(sample_name = col_names, Class = labels)
output_file <- "Rank_In/RankIn/Zm_sample_Class.txt"  
write.table(output, output_file, row.names = FALSE, sep = "\t", quote = FALSE)
