#conda activate Fisher
#conda install anaconda::pandas scipy matplotlib matplotlib-venn
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def load_and_process_data():
    """
    加载所有调控网络和显著调控网络数据，以及Chip-Seq数据，进行数据处理和交集分析。
    
    返回:
    tuple: 包含True Positives, False Positives, False Negatives, True Negatives, Odds Ratio和P-Value的元组。
    """
    all_networks = pd.read_csv('Network_GENIE3.TSV', sep='\t')
    significant_networks = pd.read_csv('Network_GENIE3_significant.txt', sep='\t')
    chip_seq_data = pd.read_csv('at_matched_results_New.txt', sep='\t')
    
    all_networks = all_networks[all_networks['TF'].isin(chip_seq_data['TF'])]   #有争议
    significant_networks = significant_networks[significant_networks['TF'].isin(chip_seq_data['TF'])]  #有争议
    filtered_chip_seq = chip_seq_data[
        chip_seq_data['target'].isin(all_networks['target']) & 
        chip_seq_data['TF'].isin(all_networks['TF'])
    ].drop_duplicates() #有争议

    # 创建调控关系集合
    all_networks_set = set(zip(all_networks['TF'], all_networks['target']))
    significant_networks_set = set(zip(significant_networks['TF'], significant_networks['target']))
    filtered_chip_seq_set = set(zip(filtered_chip_seq['TF'], filtered_chip_seq['target']))

    # 计算True Positives, False Positives, False Negatives, True Negatives
    TP = len(significant_networks_set & filtered_chip_seq_set)
    FP = len(filtered_chip_seq_set - significant_networks_set)
    FN = len(significant_networks_set - filtered_chip_seq_set)
    TN = len(all_networks_set - (significant_networks_set | filtered_chip_seq_set))

    # 构建2x2列联表
    table = [
        [TP, FP],
        [FN, TN]
    ]

    # 进行Fisher精确检验
    odds_ratio, p_value = fisher_exact(table, alternative='two-sided')
    return TP, FP, FN, TN, odds_ratio, p_value
    
def plot_venn(TP, FP, FN, p_value):
    """
    绘制基于计算得到的TP, FP, FN的Venn图，并标注P-value。
    
    参数:
    TP: int, True Positives的数量。
    FP: int, False Positives的数量。
    FN: int, False Negatives的数量。
    p_value: float, Fisher精确检验的P-value。
    """
    plt.figure(figsize=(8, 6))
    venn2(subsets=(FP, FN, TP), set_labels=('Chip-Seq Data', 'Significant Networks'))
    plt.title('Venn Diagram of Gene Regulation Networks')
    plt.figtext(0.5, 0.01, f'P-value: {p_value:.2e}', ha='center', fontsize=12)
    plt.savefig('networks_venn_diagram.pdf')
    
if __name__ == '__main__':
    """
    主程序执行区域，负责运行数据处理和分析函数，并保存相关结果。
    """
    TP, FP, FN, TN, odds_ratio, p_value = load_and_process_data()
    plot_venn(TP, FP, FN, p_value)
    with open('results.txt', 'w') as file:
        file.write(f'TP: {TP}\n')
        file.write(f'FP: {FP}\n')
        file.write(f'FN: {FN}\n')
        file.write(f'TN: {TN}\n')
        file.write(f'Odds Ratio: {odds_ratio}\n')
        file.write(f'P-Value: {p_value}\n')



# True Positives (TP): 在显著调控网络和Chip-Seq数据中都出现的调控关系。
# False Positives (FP): 在Chip-Seq数据中但不在显著调控网络中的调控关系。
# False Negatives (FN): 在显著调控网络中但不在Chip-Seq数据中的调控关系。
# True Negatives (TN): 不在显著调控网络中也不在Chip-Seq数据中的调控关系。