import pandas as pd
import numpy as np
import os
from plotnine import *
import matplotlib.pyplot as plt
import seaborn as sns
import re


def create_dfs_for_plotting(csv_downsampling, csv_coverage):
    """
    Reads and processes CSV files to prepare data for plotting.

    Args:
        csv_downsampling (str): Path to the CSV file containing downsampled sample data.
        csv_coverage (str): Path to the CSV file containing coverage information of downsampled sample data.

    Returns:
        tuple:
            - df_plot (DataFrame): Processed DataFrame excluding sex chromosomes ('chrX' and 'chrY').
            - df_cov_plot (DataFrame): Merged DataFrame with coverage information and disease labels.
            - sorted_order (list): Alphabetically sorted list of unique disease categories.

    Workflow:
    1. Reads the primary CSV file into a DataFrame `df` and removes the prefix 'Outputs/' from the 'source' column.
    2. Loads a second CSV file `df_cov` containing coverage information.
    3. Filters `df` to exclude samples from chromosomes 'chrX' and 'chrY', resulting in `df_plot`.
    4. Extracts unique 'sample' and 'disease' pairs, ensuring one entry per sample, stored in `df_dis`.
    5. Merges `df_cov` with `df_dis` on the 'sample' column to create `df_cov_plot`.
    6. Sorts unique disease categories alphabetically for consistent plotting order.
    7. Prints the number of unique samples after downsampling.

    Example:
        df_plot, df_cov_plot, sorted_order = create_df_for_plotting('data.csv')
    """
    df = pd.read_csv(csv_downsampling, index_col=0)
    df['source'] = df['source'].str.replace('Outputs/', '')
    df_cov = pd.read_csv(csv_coverage, index_col=0)

    df_plot = df[(df['chromosome'] != 'chrX') & (df['chromosome'] != 'chrY')]
    df_cov_plot = df_cov.reset_index()
    df_dis = (df[['sample', 'disease']].sort_values('sample')).groupby(['sample']).first().reset_index()
    df_cov_plot = df_cov.merge(df_dis[['sample', 'disease']], on='sample')

    # Define order for Plots
    sorted_order = sorted(df_plot['disease'].unique())  # Sort the unique categories alphabetically

    print(f"Number of total samples after downsampling: {len(df_plot['sample'].unique())}")
    return df_plot, df_cov_plot, sorted_order


def calc_chromosome_coverage(df_plot, chr_length, strand):
    """
    Calculates chromosome coverage based on median fragment data and chromosome lengths.

    Args:
        df_plot (DataFrame): DataFrame containing fragment information, including 'strand' and 'downsampling_to'.
        chr_length (str): Path to the chromosome length file in tab-separated format.
        strand (str): Strand value ('+' or '-') to filter fragments.

    Returns:
        DataFrame: DataFrame containing chromosome coverage percentages for the specified strand.
    """
    frag_chr = df_plot.loc[df_plot['strand'] == strand].groupby(by='chromosome').median('downsampling_to')

    chr_length = pd.read_csv('data/chromosome_length_diploid_genome.txt', sep='\t')
    chr_length['Length (bp)'] = chr_length['Length (bp)'].str.replace(',', '')  # Remove commas
    chr_length['Length_haploid'] = chr_length['Length (bp)'].astype(float) / 2  # Convert to float and divide by 2
    chr_length['Chromosome'] = 'chr' + chr_length['Chromosome']

    chr_length = chr_length[['Chromosome', 'Length_haploid']]
    chr_length = chr_length.set_index('Chromosome')
    frag_chr = frag_chr.merge(
        chr_length,
        how='left',
        left_index=True,
        right_index=True
    )
    frag_chr['coverage'] = frag_chr['downsampling_to'] * frag_chr['avg_frag_len'] / frag_chr['Length_haploid'] *100
    return frag_chr


"""
Visualizing Features

Function to visualize features in a Dotplot or Pointplot
    Dotplot: only unique values
    Pointplot: multiple features have the same values (frequency is shown in plot)

    feature: feature name in the column (e.g. 'CGN_NCG_80_165_ratios')
    list_df: list of dataframes containing features per chromosome and target column
"""

def visualize_features(feature, list_df):
    list_ft = []
    for df in list_df:
        cols_ft = [col for col in df.columns if col.endswith(f'{feature}')]
        # Check if we found any matching columns
        if len(cols_ft) >= 1:  # Assuming there should be at least one column
            ft = df[cols_ft]
            list_ft.append(ft)
        else:
            print(f"Warning: No {feature} columns found in dataframe")
    all_ft = pd.concat(list_ft, axis=1)
    all_ft['target'] = list_df[0]['target']
    all_ft['target'] = all_ft['target'].astype(str)
    all_ft = all_ft.reset_index()
    all_ft.columns = all_ft.columns.str.replace(f'-{feature}', '', regex=False)

    #### This block renames SampleNames
    # # Create new sample names: 'Sample1', 'Sample2', etc.
    # all_ft['SampleName'] = ['Sample' + str(i + 1) for i in range(len(all_ft))]
    # # Set the 'SampleName' back as the index

    all_ft.set_index('SampleName', inplace=True)

    # Reshape the data to long format for ggplot
    ft_long = all_ft.reset_index().melt(id_vars=['SampleName', 'target'], var_name='Chromosome', value_name=f'{feature}')

    # Calculating counts for frequency
    ft_long_counts = ft_long.groupby(['Chromosome', feature, 'target']).size().reset_index(name='count')

    # Check if all values in 'count' are 1
    if all(ft_long_counts['count'] == 1):
        # Calculate an automatic bin width based on data range
        feature_range = ft_long[f'{feature}'].max() - ft_long[f'{feature}'].min()
        num_bins = 40  # You can adjust this to make bins finer or coarser
        binwidth = feature_range / num_bins
        # Create a dot plot if all counts are 1
        plot = (
            ggplot(ft_long, aes(x='Chromosome', y=f'{feature}', fill='target')) +
            geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, binwidth=binwidth) +
            labs(title=f"Dot plot of {feature} Values across Chromosomes",
                x="Chromosomes", y=f"{feature} Values") +
            theme(axis_text_x=element_text(rotation=90, hjust=1),
                figure_size=(14, 6))
        )
    else:
        # Create a point plot with varying sizes if counts differ
        plot = (
            ggplot(ft_long_counts, aes(x='Chromosome', y=f'{feature}', color='target', size='count')) +
            geom_point(alpha=0.7) +
            labs(title=f"Point plot of {feature} Values across Chromosomes (with Frequency)",
                x="Chromosomes", y=f"{feature} Values", size="Frequency") +
            theme(axis_text_x=element_text(rotation=90, hjust=1),
                  figure_size=(12, 6))
        )

    print(plot)

"""
Visualize Fragment sizes

Function to visualize the count of all fragments sizes sequenced and aligned to reverse strand from Baseline (BL) samples
    Two options to input data:
                    1)  input_data: directory path. directory should contain the raw data csv files of reverse and BL samples
                        input_names: leave this blank. Names will automatically be the file names
                    
                    2)  input_data: list containing sample dataframes with raw features
                        input_names: list of names for the input_data as string
"""

def frag_size(input):    
    df = pd.read_csv(input, index_col=0)
    name = (input.split(sep='/')[-1]).split(sep='.')[0]
    disease = input.split(sep='/')[-2]    
    # Set chromosome as index
    df = df.set_index('chromosome')
    # Drop all columns that are not the len features
    len_pattern = r'^len_(\d+)_100k'
    cols_to_drop = [col for col in df.columns if not re.search(len_pattern, col)]
    df_len = df.drop(columns=cols_to_drop)

    # Extract the numbers in the pattern to rename columns
    new_columns = {
        col: re.search(len_pattern, col).group(1)
        for col in df_len.columns if re.search(len_pattern, col)
    }

    # Rename columns in df_len
    df_len.rename(columns=new_columns, inplace=True)

    # Group by chromosome and sum
    df_len = df_len.groupby('chromosome').sum()

    # Reshape for plotting, with 'len' and 'count' as columns
    df_long = df_len.reset_index().melt(id_vars='chromosome', var_name='len', value_name='count')
    # Convert 'len' column to float for continuous scaling
    df_long['len'] = df_long['len'].astype(float)

    # Filter for chr1 only
    df_chr1 = df_long[df_long['chromosome'] == 'chr1']

    # Plot continuous `len` values against their `count`, separated by chromosome
    plot = (
        ggplot(df_long, aes(x='len', y='count', color='chromosome'))
        + geom_line()
        + labs(x="Length", y="Count", title=f'{disease} - {name}')
        + scale_x_continuous(breaks=range(0, int(df_chr1['len'].max()) + 1, 100))
        + theme(figure_size=(15, 6))
    )
    plot_facet = plot + facet_wrap('~chromosome')  # This creates a separate panel for each chromosome
    plot_facet = plot_facet + theme(legend_position='none') # Disable the legend for color aesthetic
    plot.show()
    plot_facet.show()

    return

def sample_facet_frag_size(path, strand, chromosome):
    df_list = []  # Collect dataframes to avoid repeated merging

    # Iterate through files in the directory
    for file in os.listdir(path):
        if strand in file:
            file_path = os.path.join(path, file)
            df = pd.read_csv(file_path, index_col=0)

            # Extract name and disease from file path
            name = os.path.splitext(os.path.basename(file))[0]  # File name without extension
            disease = os.path.basename(os.path.dirname(file_path))  # Parent folder name

            # Filter len_*_100k columns
            len_pattern = r'^len_(\d+)_100k'
            len_columns = [col for col in df.columns if re.search(len_pattern, col)]
            df_len = df[len_columns].copy()

            # Rename columns to extract only numbers
            df_len.rename(
                columns={col: re.search(len_pattern, col).group(1) for col in len_columns},
                inplace=True,
            )

            # Group by chromosome, sum values, and reshape for plotting
            df_len = df_len.groupby(df['chromosome']).sum()
            df_long = (
                df_len.reset_index()
                .melt(id_vars='chromosome', var_name='len', value_name='count')
            )
            df_long['len'] = df_long['len'].astype(float)

            # Filter for the specified chromosome and format for merging
            temp_df = (
                df_long[df_long['chromosome'] == chromosome]
                .drop(columns='chromosome')
                .set_index('len')
                .rename(columns={'count': name})
            )

            df_list.append(temp_df)

    # Combine all filtered DataFrames into one
    df_chr1 = pd.concat(df_list, axis=1).reset_index()

    # Reshape for facet plot
    df_long = pd.melt(df_chr1, id_vars='len', var_name="sample", value_name="count")

    # Define the size of a single facet
    single_facet_width = 2.4  # width in inches
    single_facet_height = 1.6  # height in inches

    # Count the number of facets
    num_facets = df_long['sample'].nunique()

    # Calculate number of rows and columns for the facets
    cols =  6 # Fixed number of columns (adjust as needed)
    rows = (num_facets // cols) + (num_facets % cols > 0)  # Calculate rows based on facets and columns

    # Create the plot
    plot = (
        ggplot(df_long, aes(x='len', y='count', color='sample'))
        + geom_line()
        + facet_wrap('~sample', ncol=cols)
        + labs(x="Length", y="Count", title=f"{disease} - {chromosome}")
        + theme(figure_size=(single_facet_width * cols, single_facet_height * rows), legend_position='none')
    )
    plot.show()

# Call the function with appropriate arguments
# sample_facet_frag_size('path/to/data', 'strand_identifier', 'chr1')


def plot_frag_len(input, input_ds):
    """
    Plot Fragment Sizes before and after downsampling


    Args:
        input (str): Path to the input file containing fragment size information before downsampling.
            The file should be tab-separated with columns: ['chromosome', 'start', 'end', 'quality', 'strand'].
        input_ds (str): Path to the input file containing fragment size information after downsampling.
            The file should be tab-separated with the same columns as `input`.

    Returns:
        None: Displays two plots:
            1. KDE plot of fragment lengths before and after downsampling.
            2. Histogram of fragment lengths before and after downsampling.
    """
    size = pd.read_csv(input,
                    header=None, sep='\t')
    size_ds = pd.read_csv(input_ds,
                    header=None, sep='\t')
    size.columns = ['chromosome', 'start', 'end', 'quality', 'strand']
    size_ds.columns = ['chromosome', 'start', 'end', 'quality', 'strand']

    size['frag_length'] = size['end'] - size['start']
    size_ds['frag_length'] = size['end'] - size['start']

    chr = input.split(sep='.')[-2]

    plt.figure(figsize=(16,3))
    sns.kdeplot(data=size, x='frag_length', fill=True, color=sns.color_palette('Set2')[0])
    sns.kdeplot(data=size_ds, x='frag_length', fill=True, color=sns.color_palette('Set2')[1])
    plt.xlabel('Fragment Length')
    plt.title(f'Fragment lengths of {chr} before and after downsampling')
    plt.show()


    plt.figure(figsize=(16, 3))
    sns.histplot(size['frag_length'], bins=400, kde=True, stat="count", color=sns.color_palette('Set2')[0])  # `stat="count"` ensures y-axis is in counts
    sns.histplot(size_ds['frag_length'], bins=400, kde=True, stat="count", color=sns.color_palette('Set2')[1])
    plt.xlabel('Fragment Length')
    plt.ylabel('Count')
    plt.title(f'Fragment lengths of {chr} before and after downsampling')
    plt.show()