import pandas as pd
import os

def combine_samples(path_list):
    """
    Combines samples to one Dataframe after normalization

    Args:
        path_list (list): A list of pathes to the files to merge

    Returns:
        pd.DataFrame: One dataframe containing all samples as rows and merged chr-features as columns
    """
    df_final = None
    for path in path_list:
        base_name = os.path.basename(path)
        sample_name = base_name.split(sep='_')[0]
        strand = base_name.split(sep='_')[1]
        disease = path.split(sep='/')[-2]
        source = path.split(sep='/')[-3]
        df = pd.read_csv(path, index_col=0)
        df = df.loc[(df['avg_len'] != 0), :]
        df.reset_index(inplace=True)
        df_long = pd.melt(df, id_vars=['index'],
                        value_vars=df.columns,
                        var_name='chr-feature', value_name='norm_value')
        df_long.loc[:,'chr-feature'] = df_long['index'] + '-' + df_long['chr-feature']
        df_long = df_long[['chr-feature', 'norm_value']]
        df_wide = df_long.pivot_table(index= lambda x: sample_name, 
                        columns='chr-feature', 
                        values='norm_value', 
                        aggfunc='first')
        df_wide.index.rename('sample', inplace=True)
        df_wide['source'] = source
        # Move 'source' column to the first position
        df_wide.insert(0, 'source', df_wide.pop('source'))
        df_wide['disease'] = disease
        # Move 'disease' column to the first position
        df_wide.insert(1, 'disease', df_wide.pop('disease'))
        df_wide['strand'] = strand
        # Move 'strand' column to the second position
        df_wide.insert(2, 'strand', df_wide.pop('strand'))
        if df_final is None:
            # If df_final is not initialized, set it to the current dataframe
            df_final = df_wide.copy()
        else:
            # Incrementally append to df_final
            df_final = pd.concat([df_final, df_wide])
    return df_final