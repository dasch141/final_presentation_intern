import numpy as np
from scipy.stats import mannwhitneyu
import pandas as pd
from sklearn.model_selection import train_test_split


def move_column(df, column, position):
    """
    Moves column to a specific position in a dataframe

    Args:
        df (pd.DataFrame): DataFrame containing the column to move
        column (str): column name
        position (int): new position of column
    
    Returns:
        pd.DataFrame: DataFrame containing moved column
    """
    # Remove the target column
    target = df.pop(column)
    # Insert target as the second column (index 1)
    df.insert(position, column, target)
    return df


def del_bin_and_other_features(df, features_to_remove):
    """
    Feature cleaning - Remove non informal features

    Removed features by this function:
    * bins: bin_50_60, bin_60_70, etc. (Median, Quartile, and Mean features are enough)
    * other: list unique parts of features names to remove them
        
    Args:
        df (pd.DataFrame): The dataframe to clean
        features_to_remove (list): Unique parts of feature names

    Returns:
        pd.DataFrame: A cleaned dataframe
    """
    # Save non-informal features in a list
    bins_to_remove = [col for col in df.columns if 'bin' in col]
    # Loop trough the columns and if columnname is in the list, remove column
    df = df.drop(columns=[col for col in df.columns if any(feature in col for feature in features_to_remove)])
    df = df.drop(columns=bins_to_remove)
    return df


def selection_tpot_split(df, n=10, random_state=9):
    """
    Split samples into feature selection and tpot sets

    Args:
        df (pd.DataFrame): Dataframe containing the samples to split
        n (int): number of samples which will be selected per group
    
    Returns:
        selection_set (pd.DataFrame): contains samples for feature selection
        train_test_set (pd.DataFrame): contains remaining samples for tpot
    """
    non_cancer_samples = df[df['target'] == 0].sample(n=n, random_state=random_state)
    cancer_samples = df[df['target'] == 1].sample(n=n, random_state=random_state)
    combined_samples = pd.concat([non_cancer_samples, cancer_samples], ignore_index=False)
    list_sample_names = list(combined_samples.index)
    selection_set = df.loc[df.index.isin(list_sample_names)]
    train_test_set = df.loc[~df.index.isin(list_sample_names)]
    return selection_set, train_test_set    


def selection_tpot_split_stratified(df, feature_selection_size=0.05, random_state=9):
    """
    Split samples into feature selection and tpot sets using train_test_split from scikit learn

    Args:
        df (pd.DataFrame): contains samples for feature selection
        feature_selection_size (float): defines the size of feature selection set (e.g. 0.05 = 5% of total samples)
    
    Returns:
        pd.DataFrame: feature selection samples
        pd.DataFrame: remaining samples
    """
    X = df.drop(columns=['target']) # Features
    y  = df['target'] # Labels

    # Stratified Sampling
    X_fs, X_rest, y_fs, y_rest = train_test_split(
        X, y, test_size=(1-feature_selection_size), stratify=y, random_state=random_state
    )

    # Merge X and y back into DataFrames
    df_fs = pd.concat([X_fs, y_fs], axis=1)
    df_rest = pd.concat([X_rest, y_rest], axis=1)

    df_fs = move_column(df_fs, 'target', 2)
    df_rest = move_column(df_rest, 'target', 2)

    print(f'Feature Selection Samples: {X_fs.shape[0]}')
    print(f'Remaining Samples for Model Training: {X_rest.shape[0]}')

    return df_fs, df_rest


def mann_whitney_filter(X, y, p_value_threshold=0.05):
    """
    Perform Mann-Whitney U test on each feature to filter out non-significant ones.
    
    Parameters:
    - X: DataFrame with features.
    - y: Series or array with target values.
    - p_value_threshold: Threshold for p-values to determine significance.
    
    Returns:
    - significant_features: List of features that passed the p-value threshold.
    """
    X = X.fillna(0)
    p_values = []
    for column in X.columns:
        
        group1 = X[y == 0][column]
        group2 = X[y == 1][column]
        group1 = np.nan_to_num(group1, nan=0)
        group2 = np.nan_to_num(group2, nan=0)
        try:
            _, p_value = mannwhitneyu(group1, group2)
        except:
            p_value = 1
        p_values.append(p_value)
    
    # Create a DataFrame with features and their p-values
    p_values_df = pd.DataFrame({'feature': X.columns, 'p_value': p_values})
    
    # Filter out features with p-values greater than the threshold
    significant_features = p_values_df[p_values_df['p_value'] <= p_value_threshold]['feature'].tolist()
    return significant_features


def high_correlation_mw_filter(X_train, y_train, threshold=0.7, p_value_threshold=0.05):
    """
    Iteratively compare features and drop correlated ones based on Mann-Whitney U test.
    
    Parameters:
    - X_train: DataFrame with training features.
    - threshold: Correlation threshold to filter features.
    - p_value_threshold: Threshold for p-values to determine significance.
    
    Returns:
    - X-train: DataFrame with selected features.
    """
    significant_features = mann_whitney_filter(X_train, y_train, p_value_threshold)
    X_train = X_train[significant_features]

    remaining_features = list(X_train.columns)
    i = 0
    while i < len(remaining_features) - 1:
        print(len(remaining_features))
        print(i)
        feature1 = remaining_features[i]
        for j in range(i + 1, len(remaining_features)):
            feature2 = remaining_features[j]
            
            # Calculate the correlation coefficient
            corr_value = X_train[[feature1, feature2]].corr().iloc[0, 1]
            # Check if the features are correlated
            if abs(corr_value) > threshold:
                # Get p-values for both features
                group1_f1 = X_train[y_train == 0][feature1]
                group2_f1 = X_train[y_train == 1][feature1]
                _, p_value1 = mannwhitneyu(group1_f1, group2_f1)
                
                group1_f2 = X_train[y_train == 0][feature2]
                group2_f2 = X_train[y_train == 1][feature2]
                _, p_value2 = mannwhitneyu(group1_f2, group2_f2)
                # Compare p-values and drop the feature with higher p-value
                if p_value1 < p_value2:
                    to_drop = feature2
                elif p_value1 > p_value2:
                    to_drop = feature1
                else:  # If p-values are the same, drop randomly
                    to_drop = np.random.choice([feature1, feature2])
                
                # Remove the feature from the list and DataFrame
                remaining_features.remove(to_drop)
                X_train = X_train.drop(columns=to_drop)
                
                # Restart comparison from the most significant feature
                break
        else:
            # Only increment i if we didn't break the loop
            i += 1
    
    return X_train[remaining_features]


def feature_selection(df, threshold=0.7, p_value_threshold=0.05):
    """
    Selects features by iteratively comparing features and dropping correlated ones based on Mann-Whitney U test.

    Args:
        df (pd.DataFrame): Original DataFrame containing all features and target column.
        threshold (float): Correlation threshold to filter features.
        p_value_threshold (float): Threshold for p-values to determine significance.

    Returns:
        pd.DataFrame: DataFrame containing selected features

    """
    X = df.loc[:, ~df.columns.isin(['source', 'disease', 'target', 'strand'])]
    y = df['target']

    selected_features = high_correlation_mw_filter(X, y, threshold=threshold, p_value_threshold=p_value_threshold)
    return selected_features


def generate_export_tpot_set(df_selection, df_original, name):
    """
    Extracts selected features from the original DataFrame and exports them as a CSV file.

    Args:
        df_selection (pd.DataFrame): DataFrame containing the selected features' names as columns.
        df_original (pd.DataFrame):Original DataFrame containing all features and the target column.
        name (str): Name for the output CSV file
    
    The function creates a new DataFrame containing only the selected features and the 'target' column from the original DataFrame.
    
    It then exports this DataFrame to a CSV file located in the 'data/input_tpot/' directory with the specified file name.
    """
    feature_names = list(df_selection.columns)
    feature_names.append('target')

    input_tpot = df_original[feature_names]
    input_tpot.to_csv(f'data/input_tpot/{name}.csv')
    return