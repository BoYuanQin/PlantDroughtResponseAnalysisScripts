import networkx as nx
from itertools import combinations
def read_data(filename):
    """
    读取给定文件名的数据，返回记录的列表。
    每个记录是一个包含 'RowGene'、'ColGene'、'Identity%'、'Count_of_Gene1' 键的字典。
    """
    data = []
    with open(filename, 'r') as f:
        next(f)  # 跳过标题行
        for line in f:
            line = line.strip()
            if line:
                parts = line.split()  # 以任何空白字符分割
                if len(parts) == 4:
                    record = {
                        'RowGene': parts[0],
                        'ColGene': parts[1],
                        'Identity%': float(parts[2]),
                        'Count_of_Gene1': int(parts[3])
                    }
                    if record['RowGene'] != record['ColGene']:
                        data.append(record)
                else:
                    print(f"Warning: Line with unexpected number of columns: {line}")
    return data

def get_species(gene_id):
    """
    从基因标识符中提取物种代码。
    对应的物种返回 'Os'、'Zm' 或 'AT'。
    """
    if gene_id.startswith('Os'):
        return 'Os'
    elif gene_id.startswith('Zm'):
        return 'Zm'
    elif gene_id.startswith('AT'):
        return 'AT'
    else:
        return None  # 未知物种

def build_graph(data):
    """
    构建一个图，节点是基因，边连接同一性百分比 > 75% 的基因。
    """
    G = nx.Graph()
    for record in data:
        gene1 = record['RowGene']
        gene2 = record['ColGene']
        identity = record['Identity%']
        if identity > 75.0:
            G.add_edge(gene1, gene2, identity=identity, count=record['Count_of_Gene1'])
    return G

def find_triplets(G):
    """
    查找包含所有三个物种（Os、Zm、AT）的基因三元组，
    其中所有成对的同一性百分比均大于75%。
    返回包含三元组的列表，每个三元组包含对应的成对信息。
    """
    triplets = []
    # 获取所有大小大于等于3的团
    cliques = [clique for clique in nx.enumerate_all_cliques(G) if len(clique) >= 3]
    for clique in cliques:
        for triplet_genes in combinations(clique, 3):
            species = set(get_species(gene) for gene in triplet_genes)
            if species == {'Os', 'Zm', 'AT'}:
                gene_pairs = []
                genes = sorted(triplet_genes)  # 排序以确保顺序一致
                # 获取所有成对的基因及其属性
                for i in range(len(genes)):
                    for j in range(i+1, len(genes)):
                        gene1 = genes[i]
                        gene2 = genes[j]
                        identity = G[gene1][gene2]['identity']
                        count = G[gene1][gene2]['count']
                        gene_pairs.append({
                            'RowGene': gene1,
                            'ColGene': gene2,
                            'Identity%': identity,
                            'Count_of_Gene1': count
                        })
                triplets.append(gene_pairs)
    return triplets

def save_results(triplets, filename):
    """
    将三元组及其成对信息保存到给定的文件名，格式与原始数据相同。
    """
    with open(filename, 'w') as f:
        f.write('RowGene\tColGene\tIdentity%\tCount_of_Gene1\n')
        seen_pairs = set()
        for triplet in triplets:
            for pair in triplet:
                pair_key = tuple(sorted([pair['RowGene'], pair['ColGene']]))
                if pair_key not in seen_pairs:
                    seen_pairs.add(pair_key)
                    line = f"{pair['RowGene']}\t{pair['ColGene']}\t{pair['Identity%']}\t{pair['Count_of_Gene1']}\n"
                    f.write(line)

def TripletFinder():
    input_filename = 'readable_high_value_pairs_MCODE.txt'
    output_filename = 'triplets_results_MCODE.txt'
    data = read_data(input_filename)
    G = build_graph(data)
    triplets = find_triplets(G)
    save_results(triplets, output_filename)
    print(f"找到 {len(triplets)} 个三元组。结果已保存到 {output_filename}。")

if __name__ == '__main__':
    TripletFinder()
"""
conda install anaconda::networkx
conda activate find_gene_triplets
"""