#conda activate Three_Algorithms
import pandas as pd
import re
def extract_id_and_name(attribute):
    """
    从GFF文件的属性字段中提取基因ID和名称。
    参数:
    attribute: str, GFF文件中某个条目的属性列。
    
    返回:
    tuple: (gene_id, gene_name) 元组，包含提取的基因ID和名称。
    """
    id_match = re.search(r'ID=gene:([^;]+)', attribute)
    name_match = re.search(r'Name=([^;]+)', attribute)
    gene_id = id_match.group(1) if id_match else None
    gene_name = name_match.group(1) if name_match else None
    # print(f'Attribute: {attribute}, ID: {gene_id}, Name: {gene_name}')  # 添加调试信息
    return gene_id, gene_name

def read_gff(filename):
    """
    读取GFF文件并提取包含基因信息的条目。
    参数:
    filename: str, GFF文件的路径。
    
    返回:
    DataFrame: 包含序列名、起始位置、终止位置、基因ID和基因名称的DataFrame。
    """
    gff_data = pd.read_csv(filename, sep='\t', comment='#', header=None,
                           names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'],
                           low_memory=False)
    gff_data['seqname'] = gff_data['seqname'].astype(str)
    gff_data['start'] = pd.to_numeric(gff_data['start'], errors='coerce')
    gff_data['end'] = pd.to_numeric(gff_data['end'], errors='coerce')
    genes = gff_data.loc[gff_data['feature'] == 'gene', :].copy()
    # print(genes['attribute'].head(10))
    id_pattern = re.compile(r'ID=gene:([^;]+)')
    name_pattern = re.compile(r'Name=([^;]+)')
    # print(name_pattern)
    # genes.loc[:, 'ID'] = genes['attribute'].str.extract(id_pattern)
    # genes.loc[:, 'Name'] = genes['attribute'].str.extract(name_pattern)
    # genes[['ID', 'Name']] = genes['attribute'].apply(lambda x: pd.Series(extract_id_and_name(x)))
    genes['ID'], genes['Name'] = zip(*genes['attribute'].apply(extract_id_and_name))
    # print(genes)
    return genes[['seqname', 'start', 'end', 'ID', 'Name']]

def read_chip_seq_data(filename):
    """
    读取CHIP-seq数据文件。
    参数:
    filename: str, CHIP-seq数据文件的路径。
    
    返回:
    DataFrame: 包含染色体、起始位置、终止位置、名称、得分、链方向的DataFrame。
    """
    chip_data = pd.read_csv(filename, sep='\t', header=None,
                            names=['chr', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb'],
                            low_memory=False)
    chip_data['chr'] = chip_data['chr'].astype(str)
    chip_data['start'] = pd.to_numeric(chip_data['start'], errors='coerce')
    chip_data['end'] = pd.to_numeric(chip_data['end'], errors='coerce')
    chip_data.loc[:, 'gene_name'] = chip_data['name'].apply(lambda x: x.split(':')[0])
    return chip_data

def match_and_replace_peaks_with_gene_ids(peaks, genes):
    """
    使用基因数据更新CHIP-seq峰值数据中的基因ID。
    参数:
    peaks: DataFrame, CHIP-seq峰值数据。
    genes: DataFrame, 基因数据。
    
    返回:
    DataFrame: 更新后的CHIP-seq峰值数据。
    """
    peak_gene_ids = []
    for _, peak in peaks.iterrows():
        gene_name = peak['gene_name']
        gene_id = genes[genes['Name'] == gene_name]['ID'].values
        if len(gene_id) > 0:
            peak_gene_ids.append(gene_id[0])
        else:
            peak_gene_ids.append('Unknown')
    peaks.loc[:, 'gene_ID'] = peak_gene_ids
    return peaks

def match_peaks_to_genes(peaks, genes):
    """
    匹配CHIP-seq峰值到基因位置。
    参数:
    peaks: DataFrame, 包含峰值数据。
    genes: DataFrame, 包含基因数据。
    
    返回:
    DataFrame: 匹配结果，包含转录因子和目标基因的ID。
    """
    results = []
    for _, peak in peaks.iterrows():
        matched_genes = genes[(genes['seqname'] == peak['chr']) &
                              (genes['start'] <= peak['end']) &
                              (genes['end'] >= peak['start'])]
        if not matched_genes.empty:
            for _, gene in matched_genes.iterrows():
                results.append({
                    'TF': peak['gene_ID'],
                    'Target': gene['ID']
                })

    return pd.DataFrame(results)

def main():
    """
    主执行函数，用于整合上述步骤进行数据处理。
    """
    gff_filename = 'Arabidopsis_thaliana.TAIR10.59.gff3'
    chip_seq_filename = 'remap2022_nr_macs2_TAIR10_v1_0.bed'
    
    genes = read_gff(gff_filename)
    peaks = read_chip_seq_data(chip_seq_filename)
    peaks = match_and_replace_peaks_with_gene_ids(peaks, genes)
    matched_data = match_peaks_to_genes(peaks, genes)
    filtered_data = matched_data[matched_data['TF'] != 'Unknown']
    filtered_data.to_csv('at_matched_results.txt', sep='\t', index=False)

if __name__ == '__main__':
    main()
