import pandas as pd
from pathlib import Path
import os
import re

"""
Create directory

Function to create an directory if it does not exist yet
"""
def create_output_dir(dir):
    # Check if the directory already exist
    if not os.path.exists(dir):
        # Create directory
        os.makedirs(dir)    
    return


"""
Downsampling value

This function calculates the downsampling value needed for the 'downsampling' function
"""
def calc_ds_value(df, strand):
    # Dividing dataframe into reverse or forward strands
    df_strand = df[df['strand'].str.contains(f'{strand}')]

    # Calculate the total amount of fragments from all chromosomes per sample
    sum_frags = pd.DataFrame(df_strand.groupby('sample')['n_fragments'].sum())
    # Rename the new column to sum_fragments (total amount of fragments per sample)
    sum_frags = sum_frags.rename(columns={'n_fragments': 'sum_fragments_sample'})
    # Add the new column to the dataframe
    df_strand = df_strand.merge(sum_frags, on='sample')


    # Calculate the percentage of fragments per chromosome from total fragments in a sample
    df_strand['percent_sample'] = df_strand['n_fragments'] / df_strand['sum_fragments_sample'] * 100
    # Sort dataframe by chromosome
    df_strand = df_strand.sort_values(by=['chromosome'])

    # Determine the minimum amount of fragments between all samples per chromosome 
    min_frags = pd.DataFrame(df_strand.groupby('chromosome')['n_fragments'].min())
    # Rename the new column to downsampling_to
    min_frags = min_frags.rename(columns={'n_fragments': 'downsampling_to'})
    # Add the new column to the dataframe
    df_strand = df_strand.merge(min_frags, on='chromosome')
    # Calculate how many fragments must be removed for downsampling
    df_strand['n_remove'] = df_strand['n_fragments'] - df_strand['downsampling_to']
    # Calculate the precentage of fragments which must be removed  
    df_strand['n_remove_percent'] = df_strand['n_remove'] / df_strand['n_fragments'] * 100

    # Return the dataframe containing the downsampling informations
    return df_strand


"""
Downsampling

This function detects the minimum amount of fragments per strand for each chromosome which can be found in all samples from all diseases
and reduces all samples to the respective minimum
    e.g. afterwards the chromosome 1 in all samples of all diseases contain the same amount of fragments
Returns a Dataframe containing information for all sample chromosomes (forward and reverse strands) [frags] AND
Returns a Dataframe containing the coverage of all samples (also of removed samples) [df_coverage]
    
dir: directory name within the data directory. This directory should contain the disease directories, containing the sample directories with chromosome txt files
    e.g. 'raw_data'
coverage_threshold: set the minimum coverage which samples should have. All samples below that are removed
export_files: If False, no csv files will be exported
"""
def downsampling(dir, coverage_threshold=25, export_files=True):
    # Define input path
    input_dir = Path("data") / f'{dir}'
    # Set up columns for dataframe to store sample, chromosome, strand information and the amount of fragments
    df = pd.DataFrame({
    'sample': pd.Series(dtype='str'),
    'disease': pd.Series(dtype='str'),
    'chromosome': pd.Series(dtype='str'),
    'strand': pd.Series(dtype='str'),
    'n_fragments': pd.Series(dtype='int'),
    'avg_frag_len': pd.Series(dtype='int')
    })
    # Set up empty dictionary to store sample_name and sample data inside
    path_dict = {}
    # Set up list to store incomplete samples
    incomplete_samples = []
    # Loop through diseases
    for disease in os.listdir(input_dir):
        print(disease)
        # Loop through the samples
        disease_dir = f'{input_dir}/{disease}/'
        for sample in os.listdir(disease_dir):
            # Rename sample name if the format is not correct. Sample name should look like 'EE85759'
            sample = re.search(r'EE\d+', sample).group()
            # Print that sample was imported
            print(f'{sample} imported')
            # Loop through the files per sample
            for file_name in os.listdir(f'{disease_dir}/{sample}/'):
                # Check if sample files are still complete
                # When not, skip this sample and remove it later
                # When a file cannot be read (e.g. it is empty) it will be listed in incomplete_samples
                if sample not in incomplete_samples:
                    # Set the file path
                    file_path = f'{disease_dir}/{sample}/{file_name}'
                    # Add the name and the path to an dictionary to read the file later on
                    path_dict[file_name] = file_path
                    # Split file_name to access sample_name, disease, chromosome, and strand
                    parts = file_name.split('.')
                    try:
                        # Read in file    
                        temp_df = pd.read_csv(file_path, sep="\t", header=None)
                        # Add column names
                        temp_df.columns = ['chromosome', 'start', 'end', 'quality', 'strand']
                        # Calculate the median length of fragments
                        temp_df['avg_len'] = temp_df['end']-temp_df['start']
                        avg_len = temp_df['avg_len'].median()
                        # Fill row with sample, chromosome, strand information and the amount of fragments
                        new_row = pd.DataFrame({'sample': [parts[0]], 'disease': disease,'chromosome': [parts[2]], 'strand': [parts[-1]],
                                                'n_fragments': [temp_df.shape[0]], 'avg_frag_len': avg_len})
                        # Append row to the prepared dataframe
                        df = pd.concat([df, new_row], ignore_index=True)
                    # If file cannot be read in or information cannot be added, the sample name is added to the incomplete_samples list
                    except:
                        print(f'{sample} is incomplete')
                        # Add the sample_name to the list
                        incomplete_samples.append(parts[0])
                        # Remove duplicates
                        incomplete_samples = list(set(incomplete_samples))
    print('---')
    # Remove incomplete samples
    df = df[~df['sample'].isin(incomplete_samples)]
    # Print out which samples are removed
    for incomplete_sample in incomplete_samples:
        print(f'{incomplete_sample} was removed (incomplete)')
    # Print out when no samples are removed since all are complete
    if len(incomplete_samples) == 0:
        print('All samples imported are complete')

    # Check if coverage of sample is above the coverage_threshold
    # Calculate the median fragment length per sample
    sample_avg_len = pd.DataFrame(df.groupby('sample')['avg_frag_len'].median())
    # Calculate the total amount of fragments per sample
    sample_sum_frag = pd.DataFrame(df.groupby('sample')['n_fragments'].sum())
    # Merge the median lenght and total amount to one dataframe
    df_coverage = sample_avg_len.merge(sample_sum_frag, on='sample')
    # Set whole genome basepair number and strands analyzed
    bp_of_whole_genome = 3100000000 # whole genome has 3.1 billion (chrX and chrY are included here)
    strands = 2 # forward and reverse strand
    # Calculate the coverage per sample and add it to a new column
    df_coverage['coverage'] = df_coverage['n_fragments']*df_coverage['avg_frag_len']/bp_of_whole_genome/strands*100
    
    # Filter for samples with coverage below the threshold
    low_coverage_samples = df_coverage[df_coverage['coverage'] < coverage_threshold].index
    # Convert index containing the low coverage sample names to a list
    low_coverage_samples = low_coverage_samples.tolist()
    
    # Remove low coverage samples
    df = df[~df['sample'].isin(low_coverage_samples)]
    # Print out sample names which are removed
    for low_coverage_sample in low_coverage_samples:
        print(f'{low_coverage_sample} was removed (low coverage)')
    # Print out when all samples reached the threshold and no sample is removed
    if len(low_coverage_samples) == 0:
        print(f'All samples have a higher coverage than {coverage_threshold}')
    print('---')

    # Calculate the downsampling value per strand
    frags_rev = calc_ds_value(df, 'reverse')
    frags_fwd = calc_ds_value(df, 'forward')
    # Combine the result of both strands to one dataframe
    frags = pd.concat([frags_rev, frags_fwd], ignore_index=True)

    # Only execuded if export_files is True
    if export_files:
        # Loop trough the file_path dictionary created earlier
        for file in path_dict:
            # Split the file name to access sample_name, chromosome, and strand
            parts = file.split('.')
            # Use the file_name to access the row in the frags dataframe containing the chromosome information for that file
            selected_row = frags.loc[(frags['sample'] == parts[0]) & (frags['chromosome'] == parts[2]) & (frags['strand'] == parts[-1])]
            # Checks if information for this file was found
            if not selected_row.empty:
                # Checks if only one row of information was selected (each file should only have one row)
                if selected_row.shape[0] == 1:
                    # Extract the downsampling_value (amount of fragments to which this file is reduced)
                    downsampling_value = int(selected_row['downsampling_to'].iloc[0])
                    # Read in file
                    file_ds = pd.read_csv(path_dict[file], sep="\t", header=None)
                    # Reduce number of fragments to downsampling value (removed fragments are randomly picked)
                    file_ds = file_ds.sample(n=downsampling_value, random_state=9)
                    # Set up output_path
                    output_dir = f'{input_dir}_downsampled/{selected_row['disease'].iloc[0]}/{parts[0]}/'
                    # Create output_directory
                    create_output_dir(output_dir)
                    # Sort dataframe by its index (index numbers are still the old ones)
                    file_ds = file_ds.sort_index()
                    # Checks if downsampling worked correcly (rows in dataframe should match the downsampling_value)
                    if downsampling_value == file_ds.shape[0]:
                        # Export dataframe in the same format as input files if downsampling worked correctly
                        file_ds.to_csv(f'{output_dir}/{file}', sep='\t', header=False, index=False)
                    else:
                        print(f'Downsampling did not work correctly for {file}')
                else:
                    print('More than one downsampling value selected')
            else:
                print(f'No downsampling value exists for {file}')

    # Returns dataframe containing the chromsome/strand information (frags)
    # Returns dataframe containing the sample coverage (df_coverage)            
    return frags, df_coverage

