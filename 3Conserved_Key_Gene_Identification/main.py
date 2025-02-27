# src/main.py
import sys
import yaml
from pathlib import Path
from modules.step1_matrix import read_aligned_sequences, compute_identity_matrix, save_identity_matrix
from modules.step2_filter import CSVProcessor
from modules.step3_triplets import (
    read_data,
    build_graph,
    find_triplets,
    save_results,
    get_species
)

CONFIG_PATH = "config/setting.yaml"

def load_config():
    """加载配置文件并验证必要参数"""
    try:
        with open(CONFIG_PATH, 'r') as f:
            config = yaml.safe_load(f)
        
        parameters = config.get('parameters', {})
        required_keys = ['input_file', 'output_dir', 'output_file']
        if missing := [key for key in required_keys if key not in parameters]:
            raise ValueError(f"Missing required config keys: {', '.join(missing)}")
            
        return {
            'input_file': Path(parameters['input_file']).resolve(),
            'output_dir': Path(parameters['output_dir']).resolve(),
            'output_file': parameters['output_file'],
            'threshold': parameters.get('threshold', 75)
        }
    except FileNotFoundError:
        sys.exit(f"Config file not found: {CONFIG_PATH}")
    except Exception as e:
        sys.exit(f"Config parsing failed: {str(e)}")

def ensure_directory(path):
    """确保输出目录存在"""
    path.mkdir(parents=True, exist_ok=True)

def main():
    config = load_config()
    output_dir = config['output_dir']
    ensure_directory(output_dir)
    
    try:
        # 步骤1: 生成相似性矩阵
        matrix_output = output_dir / config['output_file']
        sequences = read_aligned_sequences(config['input_file'])
        matrix_df = compute_identity_matrix(sequences)
        save_identity_matrix(matrix_df, matrix_output)
        print(f"[1/3] 矩阵已生成: {matrix_output}")

        # 步骤2: 提取高值对
        processor = CSVProcessor(str(matrix_output))
        processor.load_csv()
        processor.extract_high_values(threshold=config['threshold'])
        
        suffix = matrix_output.stem.split('_')[-1]
        readable_output = output_dir / f"readable_high_value_pairs_{suffix}.txt"
        processor.save_readable_format(str(readable_output))
        print(f"[2/3] 高值对已提取: {readable_output}")

        # 步骤3: 查找三元组
        data = read_data(str(readable_output))
        G = build_graph(data)
        triplets = find_triplets(G)
        triplet_output = output_dir / f"triplets_results_{suffix}.txt"
        save_results(triplets, triplet_output)
        print(f"[3/3] 三元组结果已保存至: {triplet_output}")

    except Exception as e:
        sys.exit(f"处理失败: {str(e)}")

if __name__ == "__main__":
    main()
