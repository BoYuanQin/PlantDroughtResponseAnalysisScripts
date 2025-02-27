"""
conda activate Protein_Sequence_Batchs
"""
# python/find_intersection.py
import argparse
import glob
import os
from os.path import basename

def read_wgcna_key_genes(filename):
    """
    Reads the wgcna_KeyGenes.txt file and returns a set of gene names.
    """
    genes = set()
    with open(filename, 'r') as f:
        lines = f.readlines()
        # Skip the header line
        for line in lines[1:]:
            gene = line.strip()
            if gene:
                genes.add(gene)
    return genes

def read_other_key_genes(filename):
    """
    Reads other key genes files and returns a set of gene names.
    """
    genes = set()
    with open(filename, 'r') as f:
        lines = f.readlines()
        # Skip the header line
        for line in lines[1:]:
            line = line.strip()
            if line:
                columns = line.split('\t')
                gene = columns[0].strip('"')
                genes.add(gene)
    return genes

def find_intersection(genes1, genes2):
    """
    Finds the intersection of two sets of genes.
    """
    return genes1 & genes2

def save_results(filename, genes):
    """
    Saves the set of genes to a file.
    """
    with open(filename, 'w') as f:
        for gene in sorted(genes):
            f.write(f"{gene}\n")

def main():
    parser = argparse.ArgumentParser(description='Find gene intersections')
    parser.add_argument('--wgcna', required=True, help='WGCNA key genes file')
    parser.add_argument('--ppi_dir', required=True, help='PPI results directory')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    wgcna_genes = read_wgcna_key_genes(args.wgcna)
    # Get all files matching the pattern *_ppi_key_genes.txt
    files = glob.glob(os.path.join(args.ppi_dir, '*_ppi_key_genes.txt'))
    for filename in files:
        if filename == 'wgcna_KeyGenes.txt':
            continue
        other_genes = read_other_key_genes(filename)
        algorithm = basename(filename).split('_ppi_key_genes.txt')[0]
        common_genes = find_intersection(wgcna_genes, other_genes)
        # Extract algorithm name from filename
        # algorithm = filename.split('_ppi_key_genes.txt')[0]
        output_path = os.path.join(args.output_dir, f"{algorithm}_common_genes.txt")
        save_results(output_path, common_genes)
        print(f"Saved {len(common_genes)} common genes to {output_path}")

if __name__ == "__main__":
    main()
