import pandas as pd
import os
from pathlib import Path

dis_dir = '../data/raw_data_downsampled/'

for disease in os.listdir(dis_dir):
    sample_dir = os.path.join(dis_dir, disease)
    print(sample_dir)
    print(disease)
    for sample in os.listdir(sample_dir):
        print(sample)
        file_dir = os.path.join(sample_dir, sample)
        for file in os.listdir(file_dir):
            if file.endswith('ds.csv'):
                temp_file_path = os.path.join(file_dir, file)
                df = pd.read_csv(temp_file_path, index_col=0)
                df.reset_index(drop=True, inplace=True)
                df.to_csv(temp_file_path, sep='\t', index=False)
                new_file_name = '.'.join(file.split('.')[:4])
                file_path = os.path.join(file_dir, new_file_name)
                os.rename(temp_file_path, file_path)
            else:
                print(f'{file} has wrong input format')