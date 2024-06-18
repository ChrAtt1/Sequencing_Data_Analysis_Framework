from datetime import date
import numpy as np
import pandas as pd
import os
import io
import warnings
warnings.filterwarnings("ignore")


# Phase each vcf-file
def read_vcf(path):
    # Open the file
    with open(path, 'r') as f:
        # Read lines that do not start with '##' (avoiding metadata)
        lines = [l for l in f if not l.startswith('##')]

    # Parse the data into a pandas DataFrame
    return pd.read_csv(
        # Convert list of lines to a file-like object for pd.read_csv
        io.StringIO(''.join(lines)),
        # Define the data types for each column
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': float, 'FILTER': str, 'INFO': str},
        # Define the delimiter used in the file
        sep='\t'
        # Rename the '#CHROM' column to 'CHROM'
    ).rename(columns={'#CHROM': 'CHROM'})


# Phase the amplicon
def phasing_amplicon(cur_df):
    # Process phasing status and set last two columns accordingly
    for idx, row in cur_df.iterrows():
        cur_phasing_status = row['Phasing']
        ref = row['REF']
        alt = row['ALT']
        if cur_phasing_status == '0/1':
            cur_df.iat[idx, -3] = f"{ref},{alt}"
            cur_df.iat[idx, -2] = f"{ref},{alt}"
        elif cur_phasing_status == '1|0':
            cur_df.iat[idx, -3] = alt
            cur_df.iat[idx, -2] = ref
        elif cur_phasing_status == '0|1':
            cur_df.iat[idx, -3] = ref
            cur_df.iat[idx, -2] = alt
        elif cur_phasing_status in ('1|1', '1/1', '1/2'):
            cur_df.iat[idx, -3] = alt
            cur_df.iat[idx, -2] = alt

    # Set pos to row index
    cur_df = cur_df.set_index('POS', drop=False)

    # Remove columns which are not required anymore
    cur_df.drop(columns=['CHROM', 'REF', 'ALT', 'POS'], inplace=True)

    # Separate phased and unphased rows
    cur_df['Phasing temp'] = cur_df['Phasing'].str[1]
    cur_df_phased = cur_df[(cur_df['Phasing temp'] == '|') | (cur_df['Phasing temp'] == '1/1') | (cur_df['Phasing temp'] == '2/2')]
    cur_df_unphased = cur_df[(cur_df['Phasing temp'] != '|') & (cur_df['Phasing temp'] != '1/1') & (cur_df['Phasing temp'] != '2/2')]

    # Remove helper and Phasing columns
    cur_df.drop(columns=['Phasing', 'Phasing temp'], inplace=True)
    cur_df_phased.drop(columns=['Phasing', 'Phasing temp'], inplace=True)
    cur_df_unphased.drop(columns=['Phasing', 'Phasing temp'], inplace=True)

    return cur_df.T, cur_df_phased.T, cur_df_unphased.T


def ensure_unique(df, column):
    while df[column].duplicated().any():  # while there are duplicates
        # Get the positions of duplicates
        duplicate_mask = df[column].duplicated(keep='first')

        df['REF_str'] = df['REF'].astype(str)
        df['POS_str'] = df['POS'].astype(str)

        # Sort by POS and REF_str based on their lengths
        df_sorted = df.sort_values(by=['POS_str', 'REF_str'], key=lambda col: col.str.len())

        # Drop duplicates based on POS, keeping the first occurrence
        df_deduped = df_sorted.drop_duplicates(subset='POS', keep='first')

        # Drop the temporary string column
        df_deduped = df_deduped.drop(columns=['REF_str'])
        df = df_deduped.drop(columns=['POS_str'])
        df.sort_values(by=['POS'])
        df.reset_index(drop=True, inplace=True)

    return df

# Save each dataframe as excel-file
def save_dataframes_to_excel(path, prefix, suffix, dataframes):
    with pd.ExcelWriter(path + prefix + suffix + '.xlsx') as writer:
        for key, df in dataframes.items():
            df.to_excel(writer, sheet_name=key)
    #writer.close()


# Settings for the phasing tool
def setting_phasing_tool(path_vcf, path_save_result, method, amplicon_data):


    today = date.today()
    error_log = pd.DataFrame()
    vcf_files = [os.path.join(path, name) for path, _, files in os.walk(path_vcf) for name in files if
                 name.endswith('.vcf')]

    # Initialize empty dataframes for the amplicons
    cur_result = {key: pd.DataFrame() for key in amplicon_data.keys()}
    cur_results_phased = {key: pd.DataFrame() for key in amplicon_data.keys()}
    cur_results_unphased = {key: pd.DataFrame() for key in amplicon_data.keys()}
    counter_vcf = 1

    # Iterate over each vcf-file
    for cur_vcf in vcf_files:
        if(counter_vcf==36):
            x=0
        print(str(counter_vcf) + ' of ' + str(len(vcf_files)))
        counter_vcf = counter_vcf + 1
        cur_vcf = cur_vcf.replace('\\', '/')

        # Read current vcf
        cur_file = read_vcf(cur_vcf)
        # Exclude empty vcf files
        if not cur_file.empty:

            # Get DP values
            cur_file['DP'] = cur_file['INFO'].apply(lambda cur_value: int(cur_value.split('DP=')[1].split(';')[0]))

            # Get current Barcode
            cur_dataset = 'bc_' + os.path.basename(cur_vcf).split('_')[0]

            # Get amplicon-information from CHROM column
            cur_file['CHROM'] = cur_file['CHROM'].str.split('_', expand=True)[2]

            # Create new columns before dropping the old ones
            cur_file[f"{cur_dataset}_QUAL / DP"] = cur_file['QUAL'].astype(str) + ' / ' + cur_file['DP'].astype(str)
            cur_file[f"{cur_dataset}_H1"] = ""
            cur_file[f"{cur_dataset}_H2"] = ""

            cur_file['Phasing'] = cur_file.iloc[:, -5].str.split(':', expand=True)[0]

            # Remove unnecessary columns and information
            cur_file = cur_file.drop(columns=['ID', 'FILTER', 'FORMAT', 'INFO', 'QUAL', 'DP', cur_file.columns[-6]])

            # Split curFile by amplicon type
            cur_file_amplicon = {key: cur_file[cur_file['CHROM'] == key].reset_index(drop=True)
                                 for key in amplicon_data.keys()}

            # Iterate over the keys in the cur_file_amplion dictionary
            for cur_amplicon in cur_file_amplicon.keys():
                # Ensure each position is unique in the DataFrame associated with the current key
                cur_file_amplicon[cur_amplicon] = ensure_unique(cur_file_amplicon[cur_amplicon], "POS")

                # Apply the phasing_amplicon function to the current DataFrame
                cur_df, cur_df_phased, cur_df_unphased = phasing_amplicon(cur_file_amplicon[cur_amplicon])

                # If the current DataFrame in cur_result is empty, set it to the respective DataFrames
                if cur_result[cur_amplicon].empty:
                    cur_result[cur_amplicon] = cur_df
                    cur_results_phased[cur_amplicon] = cur_df_phased
                    cur_results_unphased[cur_amplicon] = cur_df_unphased
                else:
                    # Check if any of the DataFrames are empty and log error, replacing empty DataFrames with NaN values
                    for df_name, df in zip(["sample total", "sample phased", "sample unphased"],
                                           [cur_df, cur_df_phased, cur_df_unphased]):
                        if df.empty:
                            new_error = pd.DataFrame({"Error": [
                                f"Sample {cur_result[cur_amplicon].index[0].replace('_QUAL', '')}:  {df_name} of {cur_amplicon} is empty."]})
                            error_log = pd.concat([error_log, new_error], ignore_index=True)
                            df = pd.DataFrame(index=df.index)
                            df[:] = np.nan

                    # Concatenate current DataFrames to the respective result DataFrames
                    cur_result[cur_amplicon] = pd.concat([cur_result[cur_amplicon], cur_df])
                    cur_results_phased[cur_amplicon] = pd.concat([cur_results_phased[cur_amplicon], cur_df_phased])
                    cur_results_unphased[cur_amplicon] = pd.concat([cur_results_unphased[cur_amplicon], cur_df_unphased])

                    # Sort columns numerically
                    for df in [cur_result[cur_amplicon], cur_results_phased[cur_amplicon], cur_results_unphased[cur_amplicon]]:
                        df.sort_index(axis=1, inplace=True)

    # Save Error Log
    error_log = error_log.drop_duplicates(subset='Error').reset_index(drop=True)
    error_log.to_excel(path_save_result + '/' + str(today) + '_step_1_errorlog_' + method + '.xlsx', index=False)

    # Save DataFrames
    prefix = '/' + str(today) + '_step_1_vcf_import_'
    file_name = path_save_result + prefix + method + '.xlsx'
    save_dataframes_to_excel(path_save_result, prefix, method, cur_result)
    save_dataframes_to_excel(path_save_result, prefix, 'phased_' + method, cur_results_phased)
    save_dataframes_to_excel(path_save_result, prefix, 'unphased_' + method, cur_results_unphased)

    return file_name

