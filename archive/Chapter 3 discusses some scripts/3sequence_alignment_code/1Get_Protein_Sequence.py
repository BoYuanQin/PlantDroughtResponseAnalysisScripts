import pandas as pd

class FastaExtractor:
        """
    初始化FastaExtractor类的实例。
    参数:
    excel_file: str, 包含需要提取的基因标识符的Excel文件路径。
    fasta_file: str, 包含FASTA格式序列的文件路径。
    output_file: str, 用于保存提取序列的输出文件路径。
    i: int, 指定Excel文件中要读取的表格索引。
    """
    def __init__(self, excel_file, fasta_file, output_file,i):
        self.excel_file = excel_file
        self.fasta_file = fasta_file
        self.output_file = output_file
        self.identifiers = None  # 初始化时标识符设为None
        self.i = i 
    
    def read_identifiers(self):
    """
    从Excel文件读取需要提取的基因标识符。
    返回:
    identifiers: ndarray, 包含所有唯一基因标识符的数组。
    """
        # 读取Excel文件中的第一张表的第一列
        df = pd.read_excel(self.excel_file, sheet_name=self.i)  #不要求有行名，但默认有列名。
        self.identifiers = df.iloc[:, 0].unique()  # 获取第一列的所有唯一值
        return self.identifiers

    def extract_sequences(self):
    """
    根据加载的标识符从FASTA文件中提取相应的序列，并保存到输出文件。
    """
        # 如果标识符未加载，则先加载标识符
        if self.identifiers is None:
            self.read_identifiers()

        # 打开输出文件准备写入
        with open(self.output_file, 'w') as outfile:
            # 读取FASTA文件
            with open(self.fasta_file, 'r') as infile:
                write = False
                # 按行遍历FASTA文件
                for line in infile:
                    if line.startswith('>'):  # 每当遇到一个新的序列描述行
                        # 检查当前行是否包含任何标识符
                        write = any(identifier in line for identifier in self.identifiers)
                    if write:
                        # 如果当前行是一个需要写入的序列，将其写入输出文件
                        outfile.write(line)

        print("匹配的序列已经被保存到", self.output_file)

if __name__ == "__main__":
    """
    程序的主入口，用于设置和运行序列提取。
    """
    # 参数列表，每个元素是一个参数元组（excel_file, fasta_file, output_file）
    parameters = [
        ('PPI∩WGCNA基因ID-拟南芥、水稻和玉米.xlsx', 'Arabidopsis_thaliana.TAIR10.pep.all.fa', 'matched_sequences_AT.txt',0),
        ('PPI∩WGCNA基因ID-拟南芥、水稻和玉米.xlsx', 'Oryza_sativa.IRGSP-1.0.pep.all.fa', 'matched_sequences_OS.txt',1),
        ('PPI∩WGCNA基因ID-拟南芥、水稻和玉米.xlsx', 'Zea_mays.Zm-B73-REFERENCE-NAM-5.0.pep.all.fa', 'matched_sequences_ZM.txt',2)
    ]

    # 循环执行每一组参数
    for excel_file, fasta_file, output_file , i in parameters:
        extractor = FastaExtractor(excel_file, fasta_file, output_file,i)
        extractor.extract_sequences()
