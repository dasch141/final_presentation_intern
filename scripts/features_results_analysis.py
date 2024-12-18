import pandas as pd
import os
import re
from IPython.display import Markdown, display


def extract_acc_and_feat_Importance(path):
    """
    Extracting acc and weighted score values and feature importance from tpot results.

    acc: Accuracy, Sensitivity, Specificity
    
    Args:
        path (str): Path to the directory.

    Returns:
        pd.DataFrame: DataFrame containing the acc values of final classifiers (tpot thresholds were reached).
                        
                        Returns empty DataFrame if no classifiers reached tpot thresholds
    """
    df_files_acc = []
    dict_df_importance = {}
    df_acc = pd.DataFrame()
    
    for run in os.listdir(path):
        tpot_threshold_reached = False
        blind_test_found = False
        final_classifiers_path = f'{path}{run}/final_classifiers/'
        
        if not os.path.isdir(final_classifiers_path):
            print(f'{run}: final_classifiers directory does not exist')
            continue
        
        for file_name in os.listdir(final_classifiers_path):
            if 'feature_importance' in file_name:
                tpot_threshold_reached = True
                file_importance = pd.read_csv(os.path.join(final_classifiers_path, file_name))
            if 'blind_test' in file_name:
                blind_test_found = True
                file_acc = pd.read_csv(os.path.join(final_classifiers_path, file_name))
        
        # Append to df_files only if both conditions are met
        if tpot_threshold_reached and blind_test_found:
            file_acc.index = [run] * len(file_acc)
            df_files_acc.append(file_acc)

            file_importance.index = [run] * len(file_importance)
            dict_df_importance[run] = file_importance
            display(Markdown(f'**{run}: classifiers were tested on blind test set**'))
        else:
            display(Markdown((f'*{run}: tpot thresholds were not reached*')))
    
    # Only create the DataFrame if valid files exist
    if df_files_acc:
        df_acc = pd.concat(df_files_acc)
        
        df_acc.rename_axis('Feature_Set', inplace=True)
        df_acc.rename(columns={'Unnamed: 0': 'Classifier'}, inplace=True)
        df_acc = df_acc.reset_index()
        
        if 'Weighted_Score' in df_acc.columns:
            df_acc = df_acc.sort_values(by=['Weighted_Score'], ascending=False)
        else:
            print("Warning: 'Weighted_Score' column not found. Skipping sorting.")
        
    return df_acc, dict_df_importance


def filter_max_weighted_classifiers(df_acc, single_set=False):
    """
    Filters the DataFrame by keeping only rows with the maximum Weighted_Score 
    for each Feature_Set and Base_Classifier.

    Parameters:
    df_acc (pd.DataFrame): The original DataFrame containing classifiers and scores.

    Returns:
    pd.DataFrame: Filtered DataFrame with maximum Weighted_Score rows.
    """
    # Create a copy
    acc = df_acc.copy()

    # Split Classifier into two columns
    acc[['Base_Classifier', 'Classifier_ID']] = acc['Classifier'].str.rsplit('_', n=1, expand=True)

    if single_set:
        # Group by and extract rows with maximum Weighted_Score
        idx = acc.loc[acc.groupby(['Base_Classifier'])['Weighted_Score'].idxmax()].reset_index(drop=True)
        # Create the merged Classifier
        idx['Classifier'] = idx['Base_Classifier'] + '_' + idx['Classifier_ID']
        # Drop unnecessary columns
        idx = idx[['Classifier']]

        # Merge with the original DataFrame
        final_acc = pd.merge(idx, df_acc, on=['Classifier'], how='inner')
    else:
        # Group by and extract rows with maximum Weighted_Score
        idx = acc.loc[acc.groupby(['Feature_Set', 'Base_Classifier'])['Weighted_Score'].idxmax()].reset_index(drop=True)

        # Create the merged Classifier
        idx['Classifier'] = idx['Base_Classifier'] + '_' + idx['Classifier_ID']

        # Drop unnecessary columns
        idx = idx[['Feature_Set', 'Classifier']]

        # Merge with the original DataFrame
        final_acc = pd.merge(idx, df_acc, on=['Feature_Set', 'Classifier'], how='inner')

    return final_acc


def print_features_count(df):
    """
    Displays the count of non-NA values for each feature (column) in the DataFrame.

    Args:
        df (DataFrame): The DataFrame containing features to be counted.

    This function prints the name of each column and the count of its non-NA values.
    It helps to understand how many samples are present for each feature.
    """
    # Display selected features count
    print("Features selected:")
    for col in df.columns:
        count = df[col].notna().sum()
        print(f"{col}: {count}")
    return


def finding_duplicates(df):
    """
    Identifies and returns duplicate values found across different columns in the DataFrame.

    Parameters:
    - df (DataFrame): The DataFrame in which to find duplicate values.

    Returns:
    - DataFrame: A DataFrame containing:
        - 'Duplicate Value': The duplicated entry.
        - 'Row Indices': A list of row indices where the value appears.
        - 'Samples': A list of column names where the value appears.

    This function stacks the DataFrame, filters values occurring more than once,
    and groups them by the duplicate value, listing row indices and columns where they appear.
    """
    # Find all duplicates across the DataFrame
    stacked_df = df.stack().reset_index(name='Duplicate Value')

    # Filter for values appearing more than once
    duplicate_positions = (
        stacked_df[stacked_df.duplicated('Duplicate Value', keep=False)]
        .groupby('Duplicate Value')
        .agg({'level_0': list, 'level_1': list})  # Collect row indices and columns
        .reset_index()
    )

    # Rename columns for clarity
    duplicate_positions.columns = ['Duplicate Value', 'Row Indices', 'Samples']
    return duplicate_positions

def top_ten_importance_features(df):
    """
    Extracts and ranks the top ten important features from each column in a DataFrame.

    Parameters:
        df (DataFrame): A DataFrame containing feature importance scores.

    Returns:
        dfi (DataFrame): Concatenated DataFrame of sorted features for each column.
        dfiRanked (DataFrame): Intersection of the top ten features from the first two columns.
    """
    dfi = pd.DataFrame()
    list_dfs_ranked = []
    for col in df:
        df_temp = df[[col]].sort_values(col, ascending=False).reset_index()
        dfi = pd.concat([dfi, df_temp], axis=1)
        df_temp[col] = df_temp.index + 1
        df_temp = df_temp.iloc[:10]
        list_dfs_ranked.append(df_temp)
    dfiRanked = pd.merge(list_dfs_ranked[0], list_dfs_ranked[1], on ='Feature', how='inner')
    dfiRanked = dfiRanked.set_index('Feature')
    return dfi, dfiRanked

def filter_importance_features_tpot(dict_classifiers):
    """
    Filters and averages feature importance scores from multiple CSV files based on parameters count.

    Parameters:
        dict_classifiers (dict): A dictionary where keys are file paths to CSV files containing 
                                feature importance scores, and values are the corresponding 
                                parameters count to filter by.

    Returns:
        DataFrame: A concatenated DataFrame with averaged feature importance scores, 
                   renamed with classifier information.
    """
    dfi = pd.DataFrame()
    for df_path, parameters_count in dict_classifiers.items():
        classifier_name = (df_path.split(sep='_')[-3]).split(sep='-')[0]
        df = pd.read_csv(df_path, index_col=0)
        df = df.loc[df['parameters_count'] == parameters_count]
        df = pd.DataFrame(df.iloc[:, :5].mean(axis=1))
        df.rename(columns={0: f'{classifier_name}_{parameters_count}'}, inplace=True)
        dfi = pd.concat([dfi, df], axis=1)
    return dfi