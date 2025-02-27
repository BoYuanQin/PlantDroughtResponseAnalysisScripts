import pandas as pd
import scipy.stats as stats
from concurrent.futures import ProcessPoolExecutor

class DataProcessor:
    def __init__(self, file_list):
        self.file_list = file_list
    
    def process_data(self, file_name):
        data = pd.read_csv(file_name, sep='\t')
        data['Z_score'] = (data['importance'] - data['importance'].mean()) / data['importance'].std()
        data['P_value'] = stats.norm.sf(abs(data['Z_score'])) * 2
        significant_data = data[data['P_value'] < 0.05]
        # output_file = file_name.replace('.txt', '_significant.txt')
        base_name, ext = file_name.rsplit('.', 1)
        output_file = f"{base_name}_significant.txt"
        significant_data.to_csv(output_file, sep='\t', index=False)
    
    def parallel_process(self):
        with ProcessPoolExecutor() as executor:
            executor.map(self.process_data, self.file_list)
    
    def find_intersection(self):
        significant_data_list = []
        for file_name in self.file_list:
            # significant_file = file_name.replace('.txt', '_significant.txt')
            base_name, ext = file_name.rsplit('.', 1)
            significant_file = f"{base_name}_significant.txt"
            significant_data = pd.read_csv(significant_file, sep='\t')
            significant_data_list.append(significant_data)
        
        common_pairs = set(significant_data_list[0][['TF', 'target']].apply(tuple, axis=1))
        for data in significant_data_list[1:]:
            common_pairs &= set(data[['TF', 'target']].apply(tuple, axis=1))
        
        common_pairs_df = pd.DataFrame(list(common_pairs), columns=['TF', 'target'])
        merged_data = pd.concat(significant_data_list)
        intersection = pd.merge(merged_data, common_pairs_df, on=['TF', 'target'])
        
        intersection_mean = intersection.groupby(['TF', 'target']).agg({
            'importance': 'mean',
            'Z_score': 'mean',
            'P_value': 'mean'
        }).reset_index()
        output_file = f"intersection_network.txt"
        intersection_mean.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    file_list = ['network_PIDC.txt', 'Network_GRNBoost2.tsv', 'Network_GENIE3.tsv']
    processor = DataProcessor(file_list)
    processor.parallel_process()
    processor.find_intersection()
