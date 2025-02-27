import pandas as pd
import numpy as np

class CSVProcessor:
    def __init__(self, file_path):
        """Initialize the CSVProcessor with the file path."""
        self.file_path = file_path
        self.df = None
        self.final_df = None

    def load_csv(self):
        """Load a CSV file into a DataFrame."""
        # Read the CSV file, specifying tab separator and setting '--' as NaN
        self.df = pd.read_csv(
            self.file_path, 
            index_col=0
        )
        # print(self.df.head())
        
    def extract_high_values(self, threshold=80):
        """
        Extract values above a specified threshold, replacing other values with NaN.
        Also, remove rows and columns that are entirely NaN.
        """
        # Create a boolean mask where values are above the threshold
        high_values_mask = self.df > threshold
        # Replace values below the threshold with NaN
        filtered_df = self.df.where(high_values_mask)
        # Drop rows and columns that are all NaN
        self.final_df = filtered_df.dropna(how='all', axis=0).dropna(how='all', axis=1)
        print(self.final_df.head())


    def save_readable_format(self, output_path):
        """
        Save the processed DataFrame in a readable format with headers and counts of Gene1 occurrences.
        The output includes gene pairs with their identity percentage and the count of high-value pairs for each Gene1.
        """
        # Initialize list with header
        results = ["RowGene ColGene Identity% Count_of_Gene1"]

        # Create a mask for the upper triangle without the diagonal
        upper_triangle_mask = np.triu(np.ones(self.final_df.shape), k=1).astype(bool)

        # Apply the mask to get the upper triangle DataFrame
        upper_triangle = self.final_df.where(upper_triangle_mask)

        # Count the number of non-NaN values in the upper triangle for each row
        count_of_gene1 = upper_triangle.apply(lambda row: pd.notna(row).sum(), axis=1)

        # Iterate over the upper triangle and collect valid gene pairs
        for row_idx, row in enumerate(upper_triangle.index):
            for col_idx, col in enumerate(upper_triangle.columns):
                if upper_triangle_mask[row_idx, col_idx]:
                    value = upper_triangle.iat[row_idx, col_idx]
                    if pd.notna(value):
                        gene_count = count_of_gene1[row]
                        results.append(f"{row} {col} {value} {gene_count}")

        # Write the results to the output file
        with open(output_path, 'w') as file:
            for result in results:
                file.write(result + '\n')

if __name__ == "__main__":
    # Initialize the class with the CSV file path
    processor = CSVProcessor('sequence_identity_matrix_MCODE.csv')

    # Load the CSV file into a DataFrame
    processor.load_csv()

    # Extract values greater than the specified threshold and clean the data
    processor.extract_high_values(75)

    # Save the results in a readable format with headers and additional data
    processor.save_readable_format('readable_high_value_pairs_MCODE.txt')
