import pandas as pd

class CSVProcessor:
    def __init__(self, file_path):
        """
        初始化CSVProcessor类的实例。
        参数:
        file_path: str, CSV文件的路径。
        """
        self.file_path = file_path
        self.df = None
        self.final_df = None

    def load_csv(self):
         """加载CSV文件到DataFrame中。"""
        self.df = pd.read_csv(self.file_path, index_col=0)

    def extract_high_values(self, threshold=80):
    """
    提取高于指定阈值的数据子集，并移除全为NaN的行和列。
    参数:
    threshold: int, 默认为80，用于确定过滤的阈值。
    """
        high_values = self.df > threshold
        filtered_df = self.df.where(high_values)
        self.final_df = filtered_df.dropna(how='all', axis=0).dropna(how='all', axis=1)

    def save_readable_format(self, output_path):
        """Save the processed DataFrame in a readable format with headers and counts of Gene1 occurrences."""
        results = ["RowGene ColGene Identity% Count_of_Gene1"]
        count_of_gene1 = self.final_df.count(axis=1)
        for row in self.final_df.index:
            for col in self.final_df.columns:
                if pd.notna(self.final_df.at[row, col]):  
                    gene_count = count_of_gene1[row]
                    results.append(f"{row} {col} {self.final_df.at[row, col]} {gene_count}")
        """保存结果"""
        with open(output_path, 'w') as file:
            for result in results:
                file.write(result + '\n')

if __name__ == "__main__":
    """程序入口"""
    processor = CSVProcessor('924_Identity%.csv')
    processor.load_csv()
    processor.extract_high_values(70)
    processor.save_readable_format('readable_high_value_pairs.txt')
