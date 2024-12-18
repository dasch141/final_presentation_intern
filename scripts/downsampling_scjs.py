import pandas as pd
from pathlib import Path
import os
import re
from ftplib import FTP
import json
from datetime import datetime
import logging
import gc

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
Set up Logger

"""
def setup_logging(log_file):

    # Reset the root logger's handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    # Create a logger
    logger = logging.getLogger("MyLogger")
    logger.setLevel(logging.INFO) # Set the logging level

    # Create handlers
    file_handler = logging.FileHandler(log_file, mode='w')
    console_handler = logging.StreamHandler()

    # Set logging level for each handler
    file_handler.setLevel(logging.INFO)
    console_handler.setLevel(logging.INFO)

    # Create a formatter and attach it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Add the handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger

"""
Remove low chr reads samples

"""
def rmv_low_chr_reads_samples(df_strand, chr_reads_threshold):
    # Filter for samples with low reads on single chromosomes
    low_chr_reads_samples = df_strand[(df_strand['percent_sample'] < chr_reads_threshold) & (df_strand['chromosome'] != 'chrY')]['sample']
    # Convert index containing the low coverage sample names to a list
    low_chr_reads_samples = list(set(low_chr_reads_samples.tolist()))

    # Remove low chromosome reads samples
    df_strand = df_strand[~df_strand['sample'].isin(low_chr_reads_samples)]
    return df_strand, low_chr_reads_samples

"""
Downsampling value

This function calculates the downsampling value needed for the 'downsampling' function
Downsampling values are calculated by the least amount of fragments on each chromosome (strand indepentently)
"""
def calc_ds_value_chr(df, strand, chr_reads_threshold):
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


    df_strand, low_chr_reads_samples = rmv_low_chr_reads_samples(df_strand, chr_reads_threshold)

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
    return df_strand, low_chr_reads_samples


"""
Downsampling value

This function calculates the downsampling value needed for the 'downsampling' function
Downsampling values are calculated by the least amount of fragments in a sample (strand indepentently)
"""

def calc_ds_value_sample(df, strand, chr_reads_threshold):
    # Dividing dataframe into reverse or forward strands
    df_strand = df[df['strand'].str.contains(f'{strand}')]

    # Calculate the total amount of fragments from all chromosomes per sample
    sum_frags = pd.DataFrame(df_strand.groupby('sample')['n_fragments'].sum())
    # Rename the new column to sum_fragments (total amount of fragments per sample)
    sum_frags = sum_frags.rename(columns={'n_fragments': 'sum_fragments_sample'}).astype(int)
    # Add the new column to the dataframe
    df_strand = df_strand.merge(sum_frags, on='sample')


    # Calculate the percentage of fragments per chromosome from total fragments in a sample
    df_strand['percent_sample'] = (df_strand['n_fragments'] / df_strand['sum_fragments_sample'] * 100).astype('float32')
    # Sort dataframe by chromosome
    df_strand = df_strand.sort_values(by=['chromosome'])

    # Determine the minimum amount of fragments in one sample
    min_frags = pd.DataFrame(df_strand.groupby('chromosome')['sum_fragments_sample'].min())
    # Rename the new column to downsampling_to
    min_frags = min_frags.rename(columns={'sum_fragments_sample': 'downsampling_value'}).astype(int)
    # Add the new column to the dataframe
    df_strand = df_strand.merge(min_frags, on='chromosome')

    df_strand, low_chr_reads_samples = rmv_low_chr_reads_samples(df_strand, chr_reads_threshold)

    # Determine the downsampling_to value
    df_strand['downsampling_to'] = (df_strand['downsampling_value'] * df_strand['percent_sample'] / 100).round().astype(int)

    # Calculate the total amount of fragments from all chromosomes per sample
    sum_ds_frags = pd.DataFrame(df_strand.groupby('sample')['downsampling_to'].sum())
    # Rename the new column to sum_fragments (total amount of fragments per sample)
    sum_ds_frags = sum_ds_frags.rename(columns={'downsampling_to': 'sum_fragments_downsampled'}).astype(int)
    # Add the new column to the dataframe
    df_strand = df_strand.merge(sum_ds_frags, on='sample')


    # Calculate how many fragments must be removed for downsampling
    df_strand['n_remove'] = df_strand['n_fragments'] - df_strand['downsampling_to']
    # Calculate the precentage of fragments which must be removed  
    df_strand['n_remove_percent'] = (df_strand['n_remove'] / df_strand['n_fragments'] * 100).astype('float32')

    # Return the dataframe containing the downsampling informations
    return df_strand, low_chr_reads_samples



def export_downsampled_files(logger, path_dict, frags, output_downsampling):
    """
    Exports downsampled files based on the information in the frags DataFrame.

    Parameters:
        logger: Logging object for status and warnings.
        path_dict: Dictionary containing file paths to original data.
        frags: DataFrame with downsampling metadata for each file.
        output_downsampling: Directory where downsampled files will be saved.
    """
    logger.info('Start exporting downsampled files\n')
    # Loop trough the file_path dictionary created earlier
    for file_name, file_path in path_dict.items():
        # Split the file name to access sample_name, chromosome, and strand
        parts = file_name.split('.')
        # Use the file_name to access the row in the frags dataframe containing the chromosome information for that file
        selected_row = frags.loc[
            (frags['sample'] == parts[0]) &
            (frags['chromosome'] == parts[2]) &
            (frags['strand'] == parts[-1])
        ]
        # Checks if information for this file was found and only one row of information was selected (each file should only have one row)
        if not selected_row.empty and selected_row.shape[0] == 1:
            # Extract the downsampling_value (amount of fragments to which this file is reduced)
            downsampling_value = int(selected_row['downsampling_to'].iloc[0])
            # Read in file
            try:
                file_ds = pd.read_csv(file_path, sep="\t", header=None, dtype={0:object, 1:int, 2:int, 3:int, 4:object})
            except:
                file_ds = pd.read_csv(file_path, sep='\t', header=None, dtype={0:object, 1:object, 2:object, 3:object, 4:object})
                file_ds.columns = ['chromosome', 'start', 'end', 'quality', 'strand']
                file_ds = file_ds.dropna()
                file_ds = file_ds.astype({'chromosome': 'str', 'start': 'int', 'end': 'int', 'quality': 'int', 'strand': 'str'})
            # Reduce number of fragments to downsampling value (removed fragments are randomly picked)
            file_ds = file_ds.sample(n=downsampling_value, random_state=9)
            # Set up output_path
            output_dir = os.path.join(output_downsampling,selected_row['source'].iloc[0], selected_row['disease'].iloc[0], parts[0])
            # Create output_directory
            os.makedirs(output_dir, exist_ok=True)
            # Sort dataframe by its index (index numbers are still the old ones)
            file_ds = file_ds.sort_index()
            # Checks if downsampling worked correcly (rows in dataframe should match the downsampling_value)
            if downsampling_value == file_ds.shape[0]:
                # Export dataframe in the same format as input files if downsampling worked correctly
                file_ds.to_csv(os.path.join(output_dir, file_name), sep='\t', header=False, index=False)
            else:
                logger.warning(f'Downsampling did not work correctly for {file_name}')
            # Delete dataframe variable and free memory
            del file_ds
            gc.collect()

"""
Downsampling (on NAS via Mounting)

With this function, downsampling can be performed with local files or files on mounted NAS!

This function detects the minimum amount of fragments per strand for each chromosome which can be found in all samples from all diseases
and reduces all samples to the respective minimum
    e.g. afterwards the chromosome 1 in all samples of all diseases contain the same amount of fragments
Returns a Dataframe containing information for all sample chromosomes (forward and reverse strands) [frags] AND
Returns a Dataframe containing the coverage of all samples (also of removed samples) [df_coverage]
    
source_list: list of directory names within the NAS directory (nsa_dir). This directory should contain the source directories (e.g. 'Sun_2019')
nas_dir: directory path containing the source directories (mounted NAS or local)
output_downsampling: output directory for downsampled files
export_files: If False, no csv files will be exported
chr_ds_value: If True, downsampling values are calcutated due to the least amount of fragments in each chromosome
              If False, downwasmpling values are calculated due to the least amount of fragments in a sample
coverage_threshold: set the minimum coverage which samples should have. All samples below that are removed (set to None, if export_files=False)
chr_reads_threshold: threshold to exclude samples with nearly no reads on single chromosome
    (e.g. all chromosomes, except chrY, must contain at least 0.1% fragments of the total fragment amount for this sample)

"""
def downsampling_nas(source_list, ds_value, nas_dir, output_downsampling, export_files=True, coverage_threshold=25, chr_reads_threshold=0.1):
    # Set current date
    current_date = datetime.now().date().strftime("%y%m%d")
    log_file = f'results/downsampling_scjs_{current_date}_{ds_value}_log.txt'

    # Setup logging for printing logs in Terminal and log file
    logger = setup_logging(log_file)

    # Define path for local temp files (in case script crashes)
    local_temp_df = 'data/temp_files/temp_downsampling.csv'
    local_temp_file_dict = 'data/temp_files/temp_file_dict.json'

    # Set up columns for dataframe to store sample, chromosome, strand information and the amount of fragments
    df = pd.DataFrame({
    'sample': pd.Series(dtype='str'),
    'disease': pd.Series(dtype='str'),
    'source': pd.Series(dtype='str'),
    'chromosome': pd.Series(dtype='str'),
    'strand': pd.Series(dtype='str'),
    'n_fragments': pd.Series(dtype='int'),
    'avg_frag_len': pd.Series(dtype='int')
    })
    # Set up empty dictionary to store sample_name and sample data inside
    path_dict = {}
    # Set up list to store incomplete samples
    incomplete_samples = []
    for source in source_list:
        # If source data is in a subdirectory, extract the source name (e.g. Outputs/Jiang_2015 --> Jiang_2015)
        source_name = source.split(sep='/')[-1]
        logger.info(f'**{source_name}**')
        input_dir = os.path.join(nas_dir, source)
        # Loop through diseases
        for disease in os.listdir(input_dir):
            logger.info(f'*{disease}*')
            # Loop through the samples
            disease_dir = os.path.join(input_dir, disease)
            for sample in os.listdir(disease_dir):
                sample_dir=None
                if not re.search(r'\.\w{3}$', sample):
                    # Print that sample was imported
                    sample_dir = os.path.join(disease_dir, sample)
                    logger.info(f'{sample} imported')
                # Loop through the files per sample
                for file_name in os.listdir(sample_dir):
                    # Check if sample files are still complete
                    # When not, skip this sample and remove it later
                    # When a file cannot be read (e.g. it is empty) it will be listed in incomplete_samples
                    if file_name.endswith(('forward', 'reverse')) and  sample not in incomplete_samples:
                        # Split file_name to access sample_name, disease, chromosome, and strand
                        parts = file_name.split('.')
                        # Set the file path
                        file_path = os.path.join(disease_dir, parts[0], file_name)
                        # Add the name and the path to an dictionary to read the file later on
                        path_dict[file_name] = file_path

                        try:
                            # Read in file    
                            temp_df = pd.read_csv(file_path, sep="\t", header=None, dtype={0:object, 1:int, 2:int, 3:int, 4:object})                            
                            # Add column names
                            temp_df.columns = ['chromosome', 'start', 'end', 'quality', 'strand']
                            # Calculate the median length of fragments
                            temp_df['avg_len'] = temp_df['end']-temp_df['start']
                            avg_len = temp_df['avg_len'].median()
                            # Fill row with sample, chromosome, strand information and the amount of fragments
                            new_row = pd.DataFrame({'sample': [parts[0]], 'disease': disease, 'source': source_name,
                                                    'chromosome': [parts[2]], 'strand': [parts[-1]],
                                                    'n_fragments': [temp_df.shape[0]], 'avg_frag_len': avg_len})
                            # Append row to the prepared dataframe
                            df = pd.concat([df, new_row], ignore_index=True)
                        # If file cannot be read in try to remove empty rows
                        except:
                            try:
                                # Import all columns as string (e.g. float column could contain '-')
                                temp_df = pd.read_csv(file_path,sep='\t', header=None, dtype={0:object, 1:object, 2:object, 3:object, 4:object})
                                # Name columns
                                temp_df.columns = ['chromosome', 'start', 'end', 'quality', 'strand']
                                # Print rows containing empty values
                                logger.warning(f"File '{file_name}' might contain missing values. Try to remove rows containing missing values\n{temp_df[temp_df.isna().any(axis=1)]}")
                                # Drop all rows containing empty values
                                temp_df = temp_df.dropna()
                                # Change column types
                                temp_df = temp_df.astype({'chromosome': 'str', 'start': 'int', 'end': 'int', 'quality': 'int', 'strand': 'str'})
                                # Calculate the median length of fragments
                                temp_df['avg_len'] = temp_df['end']-temp_df['start']
                                avg_len = temp_df['avg_len'].median()
                                # Fill row with sample, chromosome, strand information and the amount of fragments
                                new_row = pd.DataFrame({'sample': [parts[0]], 'disease': disease, 'source': source_name,
                                                        'chromosome': [parts[2]], 'strand': [parts[-1]],
                                                        'n_fragments': [temp_df.shape[0]], 'avg_frag_len': avg_len})
                                # Append row to the prepared dataframe
                                df = pd.concat([df, new_row], ignore_index=True)
                            
                            # If file still cannot be read in or information cannot be added, the sample name is added to the incomplete_samples list
                            except:
                                logger.warning(f'{sample} is incomplete')
                                # Add the sample_name to the list
                                incomplete_samples.append(parts[0])
                                # Remove duplicates
                                incomplete_samples = list(set(incomplete_samples))
            
            # Save temp files after each finished disease
            df.to_csv(local_temp_df)
            # Save dictionary as JSON file
            with open(local_temp_file_dict, "w") as file:
                json.dump(path_dict, file)

    logger.info('All samples imported\n')
    logger.info('#####')
    # Remove incomplete samples
    df = df[~df['sample'].isin(incomplete_samples)]
    # Print out which samples are removed
    for incomplete_sample in incomplete_samples:
        logger.info(f'{incomplete_sample} was removed (incomplete)')
    # Print out when no samples are removed since all are complete
    if len(incomplete_samples) == 0:
        logger.info('All samples imported are complete')
    logger.info('---')

    # Check if coverage of sample is above the coverage_threshold

    # Set whole genome basepair number and strands analyzed
    bp_of_whole_genome = 3100000000 # whole genome has 3.1 billion (chrX and chrY are included here)
    strands = 2 # forward and reverse strand

    # Calculate the median fragment length per sample
    sample_avg_len = pd.DataFrame(df.groupby('sample')['avg_frag_len'].median())
    # Calculate the total amount of fragments per sample (for one strand)
    sample_sum_frag = pd.DataFrame((df.groupby('sample')['n_fragments'].sum())/strands)
    # Merge the median length and total amount to one dataframe
    df_coverage = sample_avg_len.merge(sample_sum_frag, on='sample')
    # Calculate the coverage per sample and add it to a new column
    df_coverage['coverage'] = df_coverage['n_fragments']*df_coverage['avg_frag_len']/bp_of_whole_genome*100

    # Filter for samples with coverage below the threshold
    low_coverage_samples = df_coverage[df_coverage['coverage'] < coverage_threshold].index
    # Convert index containing the low coverage sample names to a list
    low_coverage_samples = list(set(low_coverage_samples.tolist()))

    # Remove low coverage samples
    df = df[~df['sample'].isin(low_coverage_samples)]
    # Print out sample names which are removed
    for low_coverage_sample in low_coverage_samples:
        logger.info(f'{low_coverage_sample} was removed (low coverage)')
    # Print out when all samples reached the threshold and no sample is removed
    if len(low_coverage_samples) == 0:
        logger.info(f'All samples have coverage greater than {coverage_threshold}')
    logger.info('---')

    if ds_value == 'chr':
        logger.info('Chromosome downsampling values are calculated')
        logger.info('---')
        # Calculate the downsampling value per strand
        frags_rev, low_reads_rev = calc_ds_value_chr(df, 'reverse', chr_reads_threshold)
        frags_fwd, low_reads_fwd = calc_ds_value_chr(df, 'forward', chr_reads_threshold)
    if ds_value == 'smp':
        logger.info('Sample downsampling values are calculated')
        logger.info('---')
        # Calculate the downsampling value per strand
        frags_rev, low_reads_rev = calc_ds_value_sample(df, 'reverse', chr_reads_threshold)
        frags_fwd, low_reads_fwd = calc_ds_value_sample(df, 'forward', chr_reads_threshold)
    # Combine the result of both strands to one dataframe
    frags = pd.concat([frags_rev, frags_fwd], ignore_index=True)
    frags['percent_sample'] = frags['percent_sample'].astype('float32')
    frags['n_remove_percent'] = frags['n_remove_percent'].astype('float32')
    frags['avg_frag_len'] = frags['avg_frag_len'].astype('int32')

    # Combine the low chromosome reads lists
    low_chr_reads_samples = list(set(low_reads_rev + low_reads_fwd))
    # Remove low chromosome reads sample from coverage dataframe
    df_coverage = df_coverage[~df_coverage.index.isin(low_chr_reads_samples)]
    # Print out sample names which are removed
    for low_chr_reads_sample in low_chr_reads_samples:
        logger.info(f'{low_chr_reads_sample} was removed (low chr reads)')
    # Print out when all samples reached the threshold and no sample is removed
    if len(low_chr_reads_samples) == 0:
        logger.info(f'All samples have chr reads greater than {coverage_threshold}')
    logger.info('#####\n')

    # Save output files in results directory
    frags.to_csv(f'results/downsampling_scjs_{current_date}_{ds_value}_output.csv')
    df_coverage.to_csv(f'results/downsampling_scjs_{current_date}_{ds_value}_outputCov.csv')
    # Save dictionary containing the file paths as JSON file
    with open(f"results/downsampling_scjs_{current_date}_{ds_value}_file_dict.json", "w") as file:
        json.dump(path_dict, file)

    # Remove temp_files since script did not crash until this part
    os.remove(local_temp_df)
    os.remove(local_temp_file_dict)

    # Remove variables which are not needed for exporting
    del df, df_coverage, frags_rev, frags_fwd, local_temp_df, local_temp_file_dict, temp_df, new_row
    # Clear variables from Memory
    gc.collect()       


    # Only execuded if export_files is True
    if export_files:
        export_downsampled_files(logger, path_dict, frags, output_downsampling)
    
    logger.info('#####\n')
    try:
        # Final cleanup logs
        logger.info("All operations completed successfully.")
    finally:
        # Flush and close all log handlers
        for handler in logging.root.handlers:
            handler.flush()
            handler.close()
    

    # Returns dataframe containing the chromsome/strand information (frags)
    # Returns dataframe containing the sample coverage (df_coverage)            
    # return frags, df_coverage
















"""
Downsampling on NAS via FTP

With this function, downsampling can be performed with files on NAS using FTP!

Note: export part is still missing. Only goes through all files and calculates the downsampling value. Exports a CSV file in results directory

input_list = ['Outputs/Jiang_2015', 'Cristiano_2019', 'Snyder_2016', 'Sun_2019']
input_list = ['test_snyder']
"""
def downsampling_nas_ftp(input_list, user, password, coverage_threshold=25):
    # FTP server details
    ftp_host = "192.168.2.50"
    ftp_port = 21
    ftp_user = user
    ftp_pass = password
    base_path = '/Bioinformatics/research/FinaleDB_hg38/blood_plasma/'
    temp_path = 'data/temp_files/'


    # Set up columns for dataframe to store sample, chromosome, strand information and the amount of fragments
    df = pd.DataFrame({
    'sample': pd.Series(dtype='str'),
    'disease': pd.Series(dtype='str'),
    'source': pd.Series(dtype='str'),
    'chromosome': pd.Series(dtype='str'),
    'strand': pd.Series(dtype='str'),
    'n_fragments': pd.Series(dtype='int'),
    'avg_frag_len': pd.Series(dtype='int')
    })
    # Set up empty dictionary to store sample_name and sample data inside
    path_dict = {}
    # Set up list to store incomplete samples
    incomplete_samples = []

    # Connect to FTP to access NAS
    ftp = FTP()
    try:
        ftp.connect(ftp_host, ftp_port)
        ftp.login(ftp_user, ftp_pass)
        print("Logged in successfully!")
        
        ftp.set_pasv(True)
        ftp.cwd(base_path)
        print(ftp.pwd())    


        for input in input_list:
            disease_list = ftp.nlst(input)  # List of names in the current directory
            disease_list = [elem for elem in disease_list if not re.search(r'\.\w{3}$', elem)]
            for disease in disease_list:
                sample_list = ftp.nlst(disease)
                sample_list = [elem for elem in sample_list if not re.search(r'\.\w{3}$', elem)]
                for sample in sample_list:
                    parts = sample.split(sep='/')
                    source_id = parts[-3]
                    disease_id = parts[-2]
                    sample_id = parts[-1]
                    file_list = ftp.nlst(sample)
                    file_list = [elem for elem in file_list if not re.search(r'\.\w{3}$', elem)]
                    for file in file_list:
                        if file.endswith(('forward', 'reverse')) and sample_id not in incomplete_samples:
                            file_name = os.path.basename(file)
                            parts_file_name = file_name.split(sep='.')
                            # Add the name and the path to an dictionary to read the file later on
                            path_dict[file_name] = file
                            try:
                                # Download the file
                                remote_file = file
                                local_file = os.path.join(temp_path, file_name)
                                with open(local_file, 'wb') as file:
                                    ftp.retrbinary(f'RETR {remote_file}', file.write)
                                print(f"File '{file_name} downloaded successfully!")

                                temp_df = pd.read_csv(local_file, sep='\t', header=None)
                                # Add column names
                                temp_df.columns = ['chromosome', 'start', 'end', 'quality', 'strand']
                                # Calculate the median length of fragments
                                temp_df['avg_len'] = temp_df['end']-temp_df['start']
                                avg_len = temp_df['avg_len'].median()
                                # Fill row with sample, chromosome, strand information and the amount of fragments
                                new_row = pd.DataFrame({'sample': [sample_id], 'disease': disease_id, 'source': source_id,
                                                        'chromosome': [parts_file_name[-2]], 'strand': [parts_file_name[-1]],
                                                        'n_fragments': [temp_df.shape[0]], 'avg_frag_len': avg_len})
                                # Append row to the prepared dataframe
                                df = pd.concat([df, new_row], ignore_index=True)
                                
                                # Remove downloaded file
                                os.remove(local_file)

                            except:
                                print(f'{sample} is incomplete')
                                # Add the sample_name to the list
                                incomplete_samples.append(parts[0])
                                # Remove duplicates
                                incomplete_samples = list(set(incomplete_samples))
            
            df.to_csv(f'{temp_path}temp_result_df.csv')
            # Save to a JSON file
            with open(f"{temp_path}temp_file_dict.json", "w") as file:
                json.dump(path_dict, file)

        
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
        # Set whole genome basepair number and strands analyzed
        bp_of_whole_genome = 3100000000 # whole genome has 3.1 billion (chrX and chrY are included here)
        strands = 2 # forward and reverse strand
        # Calculate the median fragment length per sample
        sample_avg_len = pd.DataFrame(df.groupby('sample')['avg_frag_len'].median())
        # Calculate the total amount of fragments per sample
        sample_sum_frag = pd.DataFrame((df.groupby('sample')['n_fragments'].sum())/strands)
        # Merge the median lenght and total amount to one dataframe
        df_coverage = sample_avg_len.merge(sample_sum_frag, on='sample')

        # Calculate the coverage per sample and add it to a new column
        df_coverage['coverage'] = df_coverage['n_fragments']*df_coverage['avg_frag_len']/bp_of_whole_genome*100
        
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
        frags_rev = calc_ds_value_chr(df, 'reverse')
        frags_fwd = calc_ds_value_chr(df, 'forward')
        # Combine the result of both strands to one dataframe
        frags = pd.concat([frags_rev, frags_fwd], ignore_index=True)

        current_date = datetime.now().date().strftime("%y%m%d")
        frags.to_csv(f'results/downsampling_scjsFTP_output_{current_date}.csv')
        df_coverage.to_csv(f'results/downsampling_scjsFTP_outputCov_{current_date}.csv')

                    
    except Exception as e:
        print("Error:", e)
    finally:
        ftp.quit()

    """
    Upload File to NAS
    """
    # df.to_csv('data/temp_files/test.forward', sep='\t', header=None, index=False)
    # with open(local_file, 'rb') as file:
    #     ftp.storbinary('STOR output_downsampling/test.txt', file)
    # print("File uploaded successfully!")

    """
    Create directory if not exist yet
    """
    # Create directory if it not exist yet
    #os.makedirs('/path/to/directory', exist_ok=True)





