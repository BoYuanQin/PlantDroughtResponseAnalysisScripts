#conda activate Fisher
import pandas as pd
from scipy.stats import hypergeom

def calculate_hypergeometric_p_value(background_path, significant_path):
    """
    计算超几何分布的对数P值来评估显著性事件的统计显著性。
    参数:
    background_path: str, 背景集数据文件的路径，包含所有可能的事件。
    significant_path: str, 显著事件数据文件的路径，包含观察到的显著事件。
    """
    df_background = pd.read_csv(background_path, sep='\t')
    df_significant = pd.read_csv(significant_path, sep='\t')
    N = len(df_background)
    K = len(df_significant)
    # print(K )
    n = len(df_significant)
    # k = int((df_background.merge(df_significant, on=['TF', 'target']).shape[0]) / 3)
    k = 1
    log_p_value = hypergeom(N, K, n).logpmf(k)
    with open("hypergeometric_p_value.txt", "w") as file:
        file.write(f"Hypergeometric log_p_value: {log_p_value}\n")

if __name__ == "__main__":
    calculate_hypergeometric_p_value("Network_GENIE3.TSV", "Network_GENIE3_significant.txt")
