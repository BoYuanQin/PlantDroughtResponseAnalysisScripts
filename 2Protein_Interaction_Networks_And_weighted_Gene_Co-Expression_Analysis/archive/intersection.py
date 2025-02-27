# python/find_intersection.py
import argparse
import glob
import os

def main():
    parser = argparse.ArgumentParser(description='Find gene intersections')
    parser.add_argument('--wgcna', required=True, help='WGCNA key genes file')
    parser.add_argument('--ppi_dir', required=True, help='PPI results directory')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    
    wgcna_genes = read_wgcna_key_genes(args.wgcna)
    
    for ppi_file in glob.glob(os.path.join(args.ppi_dir, '*_ppi_key_genes.txt')):
        # ...处理逻辑...
        output_path = os.path.join(args.output_dir, f"{algorithm}_common_genes.txt")
        save_results(output_path, common_genes)

if __name__ == "__main__":
    main()
