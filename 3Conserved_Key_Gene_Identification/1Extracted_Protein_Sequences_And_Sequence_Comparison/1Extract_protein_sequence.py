import multiprocessing

class SequenceExtractor:
    def __init__(self, identifier_file, fasta_file, species_prefix):
        self.identifier_file = identifier_file
        self.fasta_file = fasta_file
        self.species_prefix = species_prefix
        self.identifiers = None
        self.extracted_file = f"{species_prefix}_extracted.txt"
        self.filtered_file = f"{species_prefix}_filtered.txt"

        # 获取物种对应的标识符结尾
        self.identifier_ending = self.get_identifier_ending(species_prefix)

    def get_identifier_ending(self, species_prefix):
        """根据物种前缀返回对应的标识符结尾"""
        ending_mapping = {
            'Arabidopsis_thaliana': '.1',
            'Oryza_sativa': '-01',
            'Zea_mays': 'P001',
        }
        return ending_mapping.get(species_prefix, '')

    def read_identifiers(self):
        """Reads the list of identifiers from the identifier file."""
        with open(self.identifier_file, 'r') as f:
            self.identifiers = [line.strip() for line in f]

    def extract_sequences(self):
        """Extracts sequences matching the identifiers from the FASTA file."""
        # Load identifiers if not already loaded
        if self.identifiers is None:
            self.read_identifiers()

        # Open output file to write matched sequences
        with open(self.extracted_file, 'w') as outfile:
            # Read the FASTA file
            with open(self.fasta_file, 'r') as infile:
                write = False
                # Iterate over each line in the FASTA file
                for line in infile:
                    if line.startswith('>'):  # New sequence description line
                        # Extract the gene ID from the description line
                        gene_id = None
                        for item in line.strip().split():
                            if item.startswith('gene:'):
                                gene_id = item[5:]  # Remove 'gene:' prefix
                                break
                        # Check if the gene ID is in the list of identifiers
                        if gene_id in self.identifiers:
                            write = True
                            outfile.write(line)
                        else:
                            write = False
                    elif write:
                        # Write the sequence lines corresponding to the matched gene
                        outfile.write(line)

        print(f"Matched sequences have been saved to {self.extracted_file}")

    def filter_sequences(self):
        """根据序列标识符，只保留以特定结尾的序列"""
        with open(self.extracted_file, 'r') as infile, open(self.filtered_file, 'w') as outfile:
            write_sequence = False
            for line in infile:
                if line.startswith('>'):
                    # Extract the sequence identifier
                    identifier = line[1:].split()[0].strip()
                    # Check if the identifier ends with the species-specific ending
                    if identifier.endswith(self.identifier_ending):
                        write_sequence = True
                        outfile.write(line)
                    else:
                        write_sequence = False
                elif write_sequence:
                    # Write the sequence lines corresponding to the filtered identifier
                    outfile.write(line)

        print(f"Filtered sequences have been saved to {self.filtered_file}")

def process_species(args):
    """Processes a single species by extracting and filtering sequences."""
    identifier_file, fasta_file, species_prefix = args
    extractor = SequenceExtractor(identifier_file, fasta_file, species_prefix)
    extractor.extract_sequences()
    extractor.filter_sequences()
    return extractor.filtered_file

if __name__ == "__main__":
    # List of species data: (identifier file, FASTA file, species prefix)
    species_data = [
        ('AT_fastgreedy_common_genes.txt', 'Arabidopsis_thaliana.TAIR10.pep.all.fa', 'Arabidopsis_thaliana'),
        ('OS_fastgreedy_common_genes.txt', 'Oryza_sativa.IRGSP-1.0.pep.all.fa', 'Oryza_sativa'),
        ('ZM_fastgreedy_common_genes.txt', 'Zea_mays.Zm-B73-REFERENCE-NAM-5.0.pep.all.fa', 'Zea_mays'),
    ]

    # Process all species in parallel
    with multiprocessing.Pool(processes=3) as pool:
        filtered_files = pool.map(process_species, species_data)

    # Combine the filtered sequences from all species into a single file
    combined_file = 'combined_protein_sequences.txt'
    with open(combined_file, 'w') as outfile:
        for fname in filtered_files:
            with open(fname, 'r') as infile:
                outfile.write(infile.read())

    print(f"Combined protein sequences have been saved to {combined_file}")

"""
conda activate Protein_Sequence_Batchs
"""