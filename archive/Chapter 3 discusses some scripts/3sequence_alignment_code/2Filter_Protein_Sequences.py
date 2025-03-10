# import pandas as pd

class FastaFilter:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file
    """
    从输入文件中过滤序列，只保留那些序列标识符以'1'结尾的序列。
    """
    def filter_sequences(self):
        with open(self.input_file, 'r') as infile, open(self.output_file, 'w') as outfile:
            write_sequence = False
            for line in infile:
                if line.startswith('>'):
                    # 获取第一列（序列标识符）
                    identifier = line.split('\t')[0].strip() #strip()移除空白符号
                    # 判断序列描述中的标识符是否以'.1'结尾
                    if identifier.endswith('1'):
                        write_sequence = True
                        outfile.write(line)
                    else:
                        write_sequence = False
                elif write_sequence:
                    outfile.write(line)

if __name__ == "__main__":
    """主程序入口"""
    fasta_filter = FastaFilter('matched_sequences_OS.txt', 'filtered_1_sequences_OS.txt')
    fasta_filter.filter_sequences()
