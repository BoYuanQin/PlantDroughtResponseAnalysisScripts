#!/bin/bash

# Iterate over the specified sequence files
for file in Fastgreedy_combined_protein_sequences.txt MCL_combined_protein_sequences.txt MCODE_combined_protein_sequences.txt
do
  # Check if the file exists
  if [ -f "$file" ]; then
    # Extract the algorithm name from the filename
    alg_name=$(echo "$file" | cut -d'_' -f1)
    # Define the output filename using the algorithm name
    output_file="aligned_sequences_${alg_name}.fasta"
    
    # Run MUSCLE for sequence alignment
    muscle -super5 "$file" -output "$output_file" -threads 128
    
    # Check if MUSCLE encountered an error
    if [ $? -ne 0 ]; then
      echo "Error processing $file. Moving to the next file."
      continue
    else
      echo "Successfully aligned $file. Output saved to $output_file."
    fi
  else
    # Inform the user if the file does not exist
    echo "File $file does not exist. Skipping."
  fi
done

# conda activate Sequence-Alignment