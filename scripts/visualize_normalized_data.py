import pandas as pd
import numpy as np
import os
import re
from plotnine import *
import plotly.express as px



def reshape_df_for_plotting(df, feature, color, strand):  
    """
    Visualizing Features

    Function to visualize features in a Dotplot or Pointplot
        Dotplot: only unique values
        Pointplot: multiple features have the same values (frequency is shown in plot)

        feature: feature name in the column (e.g. 'CGN_NCG_80_165_ratios')
        list_df: list of dataframes containing features per chromosome and target column
    """  
    # Select only rows containing the correct strand (forward or reverse)
    df = df.loc[df['strand'] == strand]

    # Finding columns matching the 'feature' input and adding the column names to a list
    cols_ft = [col for col in df.columns if col.endswith(f'{feature}')]
    # Check if we found any matching columns
    if len(cols_ft) >= 1:  # Assuming there should be at least one column
        # Add 'disease' to the list to also get the 'disease' information in the next step
        cols_ft.extend(['disease', 'target'])
        # Create a new dataframe containing only the selected columns listed in 'cols_ft'
        ft = df[cols_ft]
    # If no columns were found, raise a ValueError to stop the function
    else:
        raise ValueError(f"Invalid value for 'feature'. No '{feature}' columns found in dataframe")
    
    # Rename columns to only get the chromosome name (so removing the feature name)
    ft.columns = ft.columns.str.replace(f'-{feature}', '', regex=False)

    # Reshape the data to long format for plotting
    ft_long = ft.reset_index().melt(id_vars=['sample', 'disease', 'target'], var_name='chromosome', value_name=f'{feature}')
    ft_long['target'] = ft_long['target'].astype(str)
    # Group the rows by 'disease' and either 'sample' or 'chromosome'. The Median values from the non selected is used
    ft_long = pd.DataFrame(ft_long.groupby(by=['disease', color, 'target'])[feature].median()).reset_index()

    # Return the final dataframe for plotting
    return ft_long


def input_check(color, strand):
    """
    Input check

    Function to check whether color and strand input is allowed

    """
    # Validate the color input
    allowed_colors = ['chromosome', 'sample']
    if color not in allowed_colors:
        raise ValueError(f"Invalid value for 'color'. Allowed values are: {allowed_colors}")

    # Validate the color input
    allowed_strands = ['forward', 'reverse']
    if strand not in allowed_strands:
        raise ValueError(f"Invalid value for 'color'. Allowed values are: {allowed_strands}")
    return



def ggboxplot(df, feature, color, strand):
    """
    Visualizing Features

    Function to visualize features in a Boxplot using ggplot

        df: dataframe containing the chromosome-feature combination
        feature: feature name in the column (e.g. 'CGN_NCG_80_165_ratios')
        color: ['chromosome', 'sample'] defines if feature values are shown for each chromosome or each sample (median values of the other are used)
        strand: ['forward', 'reverse'] only selects files from the forward or reverse strand
    """
    input_check(color, strand)

    ft_long = reshape_df_for_plotting(df, feature, color, strand)

    plot = (
        ggplot(ft_long, aes(x='disease', y=f'{feature}')) +
        geom_boxplot(outlier_shape=None) +
        #geom_jitter(aes(color=color), height=0, width=0.35, size=2.5) +
        geom_jitter(height=0, width=0.35, size=2.5) +
        labs(title=f"Dot plot of {feature} values for each {color} across diseases",
            x="Chromosomes", y=f"{feature} Values") +
        theme(axis_text_x=element_text(),
            figure_size=(14, 4))
    )
    plot.show()

    return

def median_healthy(data, feature):
    # Calculate the median for the "Healthy" subset
    healthy_median = data[data['disease'] == 'Healthy'][feature].median()
    return healthy_median


def ggboxplot_scjs(df, feature, color, strand):

    input_check(color, strand)

    ft_long = reshape_df_for_plotting(df, feature, color, strand)
    healthy_median = median_healthy(ft_long, feature)

    # Custom labels for the legend
    custom_labels = {0: 'Non-Cancer', 1: 'Cancer'}

    plot = (
        ggplot(ft_long, aes(x='disease', y=f'{feature}')) +
        geom_boxplot(outlier_shape=None) +
        geom_jitter(aes(color='target'), height=0, width=0.2, size=1) +
        scale_color_brewer(
        type='qual',  # Qualitative palette for discrete categories
        palette='Set2',  # Choose a palette
        labels=[custom_labels[0], custom_labels[1]]  # Custom legend labels
        ) +
        # Add a vertical line using geom_segment
        geom_hline(yintercept=healthy_median, linetype='dotted', color='red') +
        labs(title=f"Dot plot of {feature} values for each {color} across diseases",
            x="Chromosomes", y=f"{feature} Values") +
        theme(axis_text_x=element_text(rotation=45, hjust=1),
            figure_size=(16, 8))
        
        
    )
    plot.show()

    return


def interactive_plot(df, feature, color, strand):

    input_check(color, strand)    

    ft_long = reshape_df_for_plotting(df, feature, color, strand)

    # Create the interactive plot
    fig = px.box(
        ft_long, 
        x='disease', 
        y=f'{feature}', 
        color=color,  # Color by chromosome
        points="all",  # Add jittered points to the boxplot
    )

    # Increase dot size for jittered points
    fig.update_traces(marker=dict(size=10))  # Apply larger dot size to jitter points

    # Customize labels
    fig.update_layout(
        xaxis_title="Diseases",
        yaxis_title=f"{feature} Values",
        title=dict(
            text=f"Interactive Dot Plot of {feature} values for each {color} across diseases",
            x=0.5,  # Center align the title
            xanchor='center',
            yanchor='top'
        ),
        margin=dict(t=40, b=20),  # Reduce top (t) and bottom (b) margins
        width=1400,  # Set the figure size
        height=400
    )

    # Calculate the median for each group
    medians = ft_long.groupby('disease')[f'{feature}'].median().reset_index()
    medians.rename(columns={f'{feature}': 'median'}, inplace=True)

    # Map categories to x-axis positions
    category_map = {category: idx for idx, category in enumerate(ft_long['disease'].unique())}

    # Add horizontal lines for medians
    for index, row in medians.iterrows():
        fig.add_shape(
            type="line",
            x0=category_map[row['disease']] - 0.4,  # Slight offset for the start of the box
            x1=category_map[row['disease']] + 0.4,  # Slight offset for the end of the box
            y0=row['median'],
            y1=row['median'],
            line=dict(color="black", width=2, dash="dot")
        )
    
    # Show the interactive plot
    fig.show()

    return


def plot_selected_features(df):
    # Assuming your DataFrame is named df
    df_long = df.reset_index().melt(
        id_vars=['sample', 'target', 'disease', 'source'], 
        var_name='feature', 
        value_name='feature_value'
    )
    df_long.rename(columns={'index': 'sample'}, inplace=True)
    df_long

    # Convert target to string for color mapping
    df_long['target'] = df_long['target'].astype(str)

    # Define custom labels for the legend
    custom_labels = ['non-cancer', 'cancer']

    # Create the plot
    plot = (
        ggplot(df_long, aes(x='feature', y='feature_value')) +
        geom_boxplot(outlier_shape=None) +
        geom_jitter(aes(color='target'), height=0, width=0.2, size=1) +
        scale_y_log10() + 
        scale_color_brewer(
            type='qual',  
            palette='Set2',  
            labels=custom_labels  
        ) +
        theme(
            axis_text_x=element_text(rotation=45, hjust=1, size=8),
            figure_size=(16, 10)
        )
    )
    plot.show()
    return


def create_plots(path_original, path_selected, path_tpot):
    df_info = pd.read_csv(path_original, index_col = 0, usecols=['sample', 'disease', 'source', 'target'])
    df_selection_set = pd.read_csv(path_selected, index_col=0)
    df_selection_set = pd.merge(df_selection_set, df_info, left_index=True, right_index=True)
    df_tpot = pd.read_csv(path_tpot, index_col=0)
    df_tpot = pd.merge(df_tpot, df_info.drop(columns=['target']), left_index=True, right_index=True)
    plot_selected_features(df_selection_set)
    plot_selected_features(df_tpot)
    return


def import_df_usecols(feature, file_path):
    df_cols = pd.read_csv(file_path, index_col=0, nrows=2)
    cols = [col for col in df_cols if col.endswith(feature)] + ['sample', 'source', 'disease', 'target', 'strand']
    df = pd.read_csv(file_path, index_col=0, usecols=cols)
    return df