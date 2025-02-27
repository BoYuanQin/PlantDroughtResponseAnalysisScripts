#!/bin/bash

# nohup iqtree -s combined_sequences_extraction2_alignment.fasta -m LG+G+F -st AA -nt AUTO -ninit 10 -ntop 5 -seed 12345 &
# conda activate Sequence-Alignment	(Includes muscle iqtree)

# Iterate over the specified sequence files
for file in AT3G15020.1_4sequences.txt AT4G26970.1_5sequences.txt
do
  # Check if the file exists
  if [ -f "$file" ]; then
    # Extract the algorithm name from the filename
    seq_name=$(echo "$file" | cut -d'_' -f1)
    # Define the output filename using the algorithm name
    output_file="aligned_sequences_${seq_name}.fasta"
    
    # Run MUSCLE for sequence alignment
    muscle -super5 "$file" -output "$output_file" -threads 128
    
	# Run iqtree for Sequence alignment Results
	iqtree -s "$output_file" -m LG+G+F -st AA -nt AUTO -bb 1000 -alrt 1000 -ninit 10 -ntop 5 -seed 12345
	
    # Check if MUSCLE encountered an error
    if [ $? -ne 0 ]; then
      echo "Error processing $file. Moving to the next file."
      continue
    else
      echo "Successfully aligned $file. Output saved to $output_file."
	  # echo "Successfully aligned $file. Output saved to $output_file."
    fi
  else
    # Inform the user if the file does not exist
    echo "File $file does not exist. Skipping."
  fi
done
