import pandas as pd
import os
from pathlib import Path

disease = 'Hepatitis_B'
input_path = f'../data/raw_data_downsampled/{disease}'
for sample in os.listdir(input_path):
    print(sample)
    sample_dir = os.path.join(input_path, sample)
    for file in os.listdir(sample_dir):
        file_path = os.path.join(sample_dir, file)
        df = pd.read_csv(file_path, sep='\t')
        if 'chromosome' in df.columns:
            df.to_csv(file_path, sep='\t', header=False, index=False)
        else:
            print(f'{file} has no column names')