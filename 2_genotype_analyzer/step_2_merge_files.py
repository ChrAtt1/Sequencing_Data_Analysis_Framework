import pandas as pd
from datetime import date

# Merge amplicon data from 2 methods
def merge_files(methods, df_gold_standard, df_file_to_compare):
    # Define constants
    INDEX_COLUMN_NAME = 'Unnamed: 0'
    method1_SUFFIX = '_' + methods[0]
    method2_SUFFIX = '_' + methods[1]

    # Nested function to sort rows based on certain column values
    def sort_rows(df):
        # Define a custom sort key function
        def custom_sort_key(index):
            parts = index.split('_')
            prefix = "_".join(parts[:-2])  # Prefix is everything except the last two elements
            suffix = "_".join(parts[-2:])  # Suffix is the last two elements

            order_dict = {
                "QUAL / DP" + method1_SUFFIX: 1,
                "H1" + method1_SUFFIX: 2,
                "H2" + method1_SUFFIX: 3,
                "QUAL / DP" + method2_SUFFIX: 4,
                "H1" + method2_SUFFIX: 5,
                "H2" + method2_SUFFIX: 6
            }

            return prefix, order_dict.get(suffix, 0)

        # Generate the sorted index with the custom sort key
        sorted_index = sorted(df.index, key=custom_sort_key)

        # Reindex the DataFrame using the sorted index
        df = df.reindex(sorted_index)

        return df

    def concat_if_column_names_match(df1, df2):
        # Extract column names from the dataframes
        df1_columns = df1.columns
        df2_columns = df2.columns

        # Extract the portion of the column name before the last "_"
        df1_columns_modified = {col.rsplit("_", 1)[0]: col for col in df1_columns}
        df2_columns_modified = {col.rsplit("_", 1)[0]: col for col in df2_columns}

        # Find the intersection of the keys
        matching_keys = df1_columns_modified.keys() & df2_columns_modified.keys()

        # Get the actual column names from the keys
        matching_columns_df1 = [df1_columns_modified[key] for key in matching_keys]
        matching_columns_df2 = [df2_columns_modified[key] for key in matching_keys]

        # Select only the columns that match and concatenate the dataframes
        df_concat = pd.concat([df1[matching_columns_df1], df2[matching_columns_df2]], axis=0)
        return df_concat

    # Function to join non-null values
    def join_non_nulls(series):
        return ', '.join(series.dropna().astype(str))

    # Add suffixes to INDEX_COLUMN_NAME in the respective dataframes
    df_gold_standard[INDEX_COLUMN_NAME] += method1_SUFFIX
    df_file_to_compare[INDEX_COLUMN_NAME] += method2_SUFFIX

    # Set INDEX_COLUMN_NAME as index for both dataframes
    df_gold_standard.set_index(INDEX_COLUMN_NAME, inplace=True)
    df_file_to_compare.set_index(INDEX_COLUMN_NAME, inplace=True)

    # Transpose both dataframes
    df_gold_standard_T = df_gold_standard.transpose()
    df_file_to_compare_T = df_file_to_compare.transpose()

    # Concatenate the transposed dataframes
    df_concatenated_T = concat_if_column_names_match(df_gold_standard_T, df_file_to_compare_T)

    # Group by index (which corresponds to the row names in the original dataframes) and aggregate using the custom function
    df_concatenated_T = df_concatenated_T.groupby(level=0).agg(join_non_nulls)

    # Transpose the dataframe back to original form
    df_concatenated = df_concatenated_T.transpose()
    df_concatenated.columns = df_concatenated.columns.astype(int)

    df_concatenated_sorted = sort_rows(df_concatenated)

    # Rename the Unnamed: 0 column
    df_concatenated_sorted = df_concatenated_sorted.rename(columns={INDEX_COLUMN_NAME: 'variant position by provided reference sequence'})

    return df_concatenated_sorted


# setting for merging tool
def merging_tool_setting(amplicon_data, methods, path_file_method1, path_file_method2, path_output_excel):
    # Get the current date
    today = date.today()

    # Define the name of the output Excel file
    file_name_merged = path_output_excel + '/' + str(today) + '_step_2_' + methods[0] + '_and_' + methods[1] + '_merged.xlsx'

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer_excel = pd.ExcelWriter(file_name_merged, engine='xlsxwriter')

    # Iterate through the keys in amplicon_data to perform the merging
    for sheet_name in amplicon_data.keys():
        print('Merging: ' + sheet_name)

        # Load data from the current sheet
        df_method1 = pd.read_excel(path_file_method1, sheet_name=sheet_name)

        # Load data from the current sheet
        df_method2 = pd.read_excel(path_file_method2, sheet_name=sheet_name)

        # Merge df_gold_standard with df_to_compare and update df_gold_standard
        df_merged = merge_files(methods, df_method1, df_method2)

        # Write df_gold_standard_export to the current sheet in the output Excel file
        df_merged.to_excel(writer_excel, sheet_name=sheet_name)

    # Close the Pandas Excel writer
    writer_excel.close()

    return file_name_merged
