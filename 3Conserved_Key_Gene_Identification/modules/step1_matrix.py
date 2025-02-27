import pandas as pd
import numpy as np
from Bio import SeqIO

def read_aligned_sequences(filename):
    """
    Read aligned sequences from a FASTA file.
    Return a dictionary with sequence IDs as keys and sequences as values.
    """
    sequences = {}
    try:
        with open(filename, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                sequences[record.id] = str(record.seq)
        if not sequences:
            print("No sequences found in the file.")
    except Exception as e:
        print(f"Error reading file: {e}")
    return sequences

def calculate_identity(seq1, seq2):
    """
    Calculate the identity percentage between two aligned sequences, excluding gaps.
    """
    matches = 0
    length = 0
    for a, b in zip(seq1, seq2):    #长度要一致
        # Exclude positions with gaps in either sequence
        if a != '-' and b != '-':
            length += 1
            if a == b:
                matches += 1
    if length == 0:
        identity = 0
    else:
        identity = (matches / length) * 100
    return identity

def compute_identity_matrix(sequences):
    """
    Compute the identity matrix for all sequences.
    Return a pandas DataFrame.
    """
    ids = list(sequences.keys())
    n = len(ids)
    identity_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            seq1 = sequences[ids[i]]
            seq2 = sequences[ids[j]]
            identity = calculate_identity(seq1, seq2)
            identity_matrix[i][j] = identity
            identity_matrix[j][i] = identity  # Symmetric matrix
    df = pd.DataFrame(identity_matrix, index=ids, columns=ids)
    return df

def save_identity_matrix(df, output_filename):
    """
    Save the identity matrix to a CSV file.
    """
    try:
        df.to_csv(output_filename)
        print(f"Identity matrix has been saved to {output_filename}")
    except Exception as e:
        print(f"Error saving file: {e}")

if __name__ == "__main__":
    # Input and output filenames
    input_filename = 'aligned_sequences_MCODE.fasta'
    output_filename = 'sequence_identity_matrix_MCODE.csv'

    # Read aligned sequences from the file
    sequences = read_aligned_sequences(input_filename)

    # Check if any sequences were read
    if sequences:
        # Compute the identity matrix
        identity_df = compute_identity_matrix(sequences)

        # Save the identity matrix to a CSV file
        save_identity_matrix(identity_df, output_filename)
    else:
        print("No sequences to process.")
"""
conda activate Protein_Sequence_Batchs
"""