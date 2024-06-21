import math
from collections import Counter
from datetime import date
import scipy as sp
import numpy as np
import pandas as pd
import re


# Check if value is nan.
def check_is_nan(import_value):
    try:
        return math.isnan(float(import_value))
    except:
        return False


# Cut off columns which are outside of the limits.
def cut_off_sequence(df_to_analyze, cutoff_start, cutoff_end):
    # Separate the first column and the remaining columns
    df_to_analyze_temp = df_to_analyze.iloc[:, :1]
    df_to_analyze = df_to_analyze.iloc[:, 1:]

    # Filter the columns based on the cutoff start and end positions
    df_to_analyze = df_to_analyze.loc[:, (df_to_analyze.columns.astype(float) > float(cutoff_start)) & (
            df_to_analyze.columns.astype(float) < float(cutoff_end))]

    # Concatenate the first column with the filtered DataFrame
    df_to_analyze = pd.concat([df_to_analyze_temp, df_to_analyze], axis=1)

    return df_to_analyze


# Set values to nan if dp or qual is under or over the limit
def filter_data_qual_dp(df_to_analyze, methods, cutoff_qual_method1, cutoff_dp_method1, cutoff_qual_method2,
                        cutoff_dp_method2):
    print('## Filter Data Qual and DP ##')

    # Extract method1 and method2
    method1 = methods[0]
    method2 = methods[1]

    # Iterates though every column and row to set haplotype to nan if they are under the cutoff sequence
    for columnIdx, columnData in enumerate(df_to_analyze.items()):
        for idxRow, curRow in enumerate(columnData[1]):
            # check if Qual / DP is in Name and that the value is not nan (= 3 letters)
            if 'QUAL / DP_' + method2 in df_to_analyze.index[idxRow] and len(str(curRow)) > 3:
                # extract cur Qual and DP
                cur_qual = float(curRow.split('/')[0])
                cur_dp = float(curRow.split('/')[1])

                # if under condition set QUAL / DP, H1 and H2 to nan
                if cur_qual < float(cutoff_qual_method2) or cur_dp < float(cutoff_dp_method2):
                    df_to_analyze.loc[:, columnData[0]][idxRow] = np.nan
                    df_to_analyze.loc[:, columnData[0]][idxRow + 1] = np.nan
                    df_to_analyze.loc[:, columnData[0]][idxRow + 2] = np.nan

            # check if Qual / DP for method1 is in Name and that the value is not nan (= 3 letters)
            elif 'QUAL / DP_' + method1 in df_to_analyze.index[idxRow] and len(str(curRow)) > 3:
                # extract current Qual and DP for method1
                cur_qual_method1 = float(curRow.split('/')[0])
                cur_dp_method1 = float(curRow.split('/')[1])

                # if under condition set QUAL / DP, H1 and H2 to nan for method1
                if cur_qual_method1 < float(cutoff_qual_method1) or cur_dp_method1 < float(cutoff_dp_method1):
                    df_to_analyze.loc[:, columnData[0]][idxRow] = np.nan
                    df_to_analyze.loc[:, columnData[0]][idxRow + 1] = np.nan
                    df_to_analyze.loc[:, columnData[0]][idxRow + 2] = np.nan

    # Return the filtered and cutting Dataframe
    return df_to_analyze


# Generate combinations from given dataframe and positions.
def get_combinations(df_to_analyze, methods, positions, cur_amplicon):
    # Extract method names
    method1 = methods[0]
    method2 = methods[1]

    # Prepare a list to store the existing combinations
    existing_combinations = []

    # Special handling based on the 'Sequencing Direction', reverse the positions if it's 'backward'
    # If the sequencing direction is 'backward', the positions list is reversed.
    if cur_amplicon['Sequencing Direction'] == 'backward':
        positions = positions[::-1]

    # Iterate over each row in the dataframe
    for idx_row, cur_row in df_to_analyze.iterrows():
        # If the row index contains '_H'
        if '_H' in idx_row:
            # Determine the method based on the row index
            if method1 in idx_row:
                method = method1
            elif method2 in idx_row:
                method = method2

            # Generate combination from non-nan values
            # For each position and corresponding entry in the row, if the entry is not nan,
            # create a string of the format "positionEntry " and join all such strings.
            cur_haplo_combination = ''.join([f"{int(position)}{curEntry} " for position, curEntry
                                             in zip(positions, cur_row) if not check_is_nan(curEntry)])

            # Append the combination to list along with its method
            existing_combinations.append((cur_haplo_combination, method))

    # Filter out empty combinations
    existing_combinations = list(filter(None, existing_combinations))

    # Count the occurrence of each combination
    # This is a dictionary where the keys are combinations and the values are their counts.
    count_combinations = Counter(existing_combinations)

    # Sort the combinations by their counts in descending order and convert to a list
    list_output = count_combinations.most_common()

    # Return the sorted list of combinations and their counts
    return list_output

# Generate combinations from given dataframe and positions.
def extract_haplotype_frequencies(df_to_analyze, methods):
    print('## Extract Halpotype Frequencies ##')

    # Extract methods
    method1 = methods[0]
    method2 = methods[1]

    # Check and summarize list
    def check_and_sum(lst):
        return sum(item for sublist in lst for item in sublist if isinstance(item, (int, float)))

    # Check if value is empty
    def check_is_empty(import_value, idx):
        try:
            return import_value[idx][0]
        except:
            return ''

    def add_if_not_empty(func_return):
        return '/' + str(func_return) if func_return != '' else ''

    def calculate_frequency(lst1, lst2):
        # Concatenate the two lists
        lst = lst1 + lst2

        # Initialize a Counter to keep track of the character counts
        counter = Counter()

        # Total count of all characters
        total_count = 0

        # Iterate over the list of tuples
        for tup in lst:
            # The string of characters is the first item in the tuple
            chars = tup[0]

            # The count is the second item in the tuple
            count = tup[1]

            # Multiply the count by the number of occurrences of each character in the string
            for char in chars:
                counter[char] += count

            # Add the count to the total count
            total_count += count * len(chars)

        # Calculate the frequency of each character
        frequency = {char: round(count / total_count, 4) for char, count in counter.items()}

        return frequency

    def check_entry(entry):
        # If the value of the entry is None, create a new entry with key "None" and value 0
        if entry is None:
            entry = ("None", 0)
        return entry

    # Check if value is number
    def check_is_number(import_value, idx_1, idx_2):
        try:
            return import_value[idx_1][idx_2]
        except:
            return 0

    def rounded_expected_values(expected_values, target_sum):
        rounded_values = [int(round(val)) for val in expected_values]
        correction = target_sum - sum(rounded_values)
        rounded_values[-1] += correction
        return rounded_values

    # Make a copy of the input DataFrame and remove the last three columns
    df_to_analyze_return = df_to_analyze
    df_to_analyze = df_to_analyze.iloc[:, :-6]
    row_names = list(df_to_analyze.index)

    # hardy weinberg is calculated
    # no_of_rows = 31
    # hardy_weinberg = 1

    # hardy weinberg is not calculated
    no_of_rows = 18
    hardy_weinberg = 0

    df_to_analyze_return = insert_empty_rows_method1_method2(df_to_analyze_return, methods, no_of_rows, hardy_weinberg)

    # Replace NaN values in the DataFrame with the value from the first row of the same column
    for idxCol, curColumn in enumerate(df_to_analyze.items()):
        for idxRow, curRow in enumerate(curColumn[1]):
            if 'H' in str(row_names[idxRow]) and check_is_nan(curRow):
                df_to_analyze.iloc[idxRow, idxCol] = df_to_analyze.iloc[2, idxCol]

    # Iterate through each column in the DataFrame
    for idxCol, curColumn in enumerate(df_to_analyze.items()):

        # Initialize lists to store homozygous and heterozygous haplotype frequencies for method1 and method2
        cur_homozygot_method1, cur_heterozygot_method1, cur_homozygot_method2, cur_heterozygot_method2 = [], [], [], []

        # Iterate through each row in the current column
        for idxRow, curRow in enumerate(curColumn[1]):
            # Process method1 homozygous and heterozygous haplotype frequencies
            if 'H1' in str(row_names[idxRow]) and method1 in str(row_names[idxRow]) and len(str(curRow)) <= 3:
                # Homozygous case
                if curRow == curColumn[1].iloc[idxRow + 1] and len(str(curRow)) <= 1:
                    cur_homozygot_method1.append(''.join(sorted(curRow + curColumn[1].iloc[idxRow + 1])))
                # Heterozygous case
                elif curRow == curColumn[1].iloc[idxRow + 1] and len(str(curRow)) > 1:
                    cur_heterozygot_method1.append(''.join(sorted(set(curRow.replace("/", "")))))
                else:
                    cur_heterozygot_method1.append(''.join(sorted(curRow + curColumn[1].iloc[idxRow + 1])))

            # Process method2 homozygous and heterozygous haplotype frequencies
            elif 'H1' in str(row_names[idxRow]) and method2 in str(row_names[idxRow]) and len(str(curRow)) <= 3:
                # Homozygous case
                if curRow == curColumn[1].iloc[idxRow + 1] and len(str(curRow)) <= 1:
                    cur_homozygot_method2.append(''.join(sorted(curRow + curColumn[1].iloc[idxRow + 1])))
                # Heterozygous case
                elif curRow == curColumn[1].iloc[idxRow + 1] and len(str(curRow)) > 1:
                    cur_heterozygot_method2.append(''.join(sorted(set(curRow.replace("/", "")))))
                else:
                    cur_heterozygot_method2.append(''.join(sorted(curRow + curColumn[1].iloc[idxRow + 1])))

        # Count the occurrences of each haplotype combination for both method1 and method2 datasets
        list_output_cur_homozygot_method1 = Counter(cur_homozygot_method1).most_common()
        list_output_cur_heterozygot_method1 = Counter(cur_heterozygot_method1).most_common()
        list_output_cur_homozygot_method2 = Counter(cur_homozygot_method2).most_common()
        list_output_cur_heterozygot_method2 = Counter(cur_heterozygot_method2).most_common()

        # Calculate the total haplotype frequencies for method1 and method2 datasets
        total_method1 = check_and_sum(list_output_cur_homozygot_method1) + check_and_sum(
            list_output_cur_heterozygot_method1)
        total_method2 = check_and_sum(list_output_cur_homozygot_method2) + check_and_sum(
            list_output_cur_heterozygot_method2)

        homozygot1_observed_method1 = check_is_number(list_output_cur_homozygot_method1, 0, 1)
        heterozygot_observed_method1 = check_is_number(list_output_cur_heterozygot_method1, 0, 1)
        homozygot2_observed_method1 = check_is_number(list_output_cur_homozygot_method1, 1, 1)
        homozygot_alt_method1 = check_is_number(list_output_cur_homozygot_method1, 2, 1) \
                                + check_is_number(list_output_cur_homozygot_method1, 3, 1)
        heterozygot_alt_method1 = check_is_number(list_output_cur_heterozygot_method1, 1, 1) \
                                  + check_is_number(list_output_cur_heterozygot_method1, 2, 1)

        # Insert the haplotype frequency values into the return DataFrame
        df_to_analyze_return.iloc[len(df_to_analyze) + 1, idxCol] = homozygot1_observed_method1
        df_to_analyze_return.iloc[len(df_to_analyze) + 2, idxCol] = heterozygot_observed_method1
        df_to_analyze_return.iloc[len(df_to_analyze) + 3, idxCol] = homozygot2_observed_method1
        df_to_analyze_return.iloc[len(df_to_analyze) + 4, idxCol] = homozygot_alt_method1
        df_to_analyze_return.iloc[len(df_to_analyze) + 5, idxCol] = heterozygot_alt_method1

        curCategory_method1 = ''
        curCategory_method1 += add_if_not_empty(check_is_empty(list_output_cur_homozygot_method1, 0))
        curCategory_method1 += add_if_not_empty(check_is_empty(list_output_cur_heterozygot_method1, 0))
        curCategory_method1 += add_if_not_empty(check_is_empty(list_output_cur_homozygot_method1, 1))

        # Remove the leading backslash
        if curCategory_method1 and curCategory_method1[0] == '/':
            curCategory_method1 = curCategory_method1[1:]

        frequency_method1 = calculate_frequency(list_output_cur_homozygot_method1,
                                                list_output_cur_heterozygot_method1)

        items_method1 = iter(frequency_method1.items())
        first_item_method1 = next(items_method1, None)
        second_item_method1 = next(items_method1, None)
        third_item_method1 = next(items_method1, None)
        first_item_method1 = check_entry(first_item_method1)
        second_item_method1 = check_entry(second_item_method1)
        third_item_method1 = check_entry(third_item_method1)
        first_item_method1 = first_item_method1[1]
        second_item_method1 = second_item_method1[1] + third_item_method1[1]

        df_to_analyze_return.iloc[len(df_to_analyze) + 6, idxCol] = first_item_method1
        df_to_analyze_return.iloc[len(df_to_analyze) + 7, idxCol] = second_item_method1
        df_to_analyze_return.iloc[len(df_to_analyze) + 8, idxCol] = curCategory_method1

        # method 2
        homozygot1_observed_method2 = check_is_number(list_output_cur_homozygot_method2, 0, 1)
        heterozygot_observed_method2 = check_is_number(list_output_cur_heterozygot_method2, 0, 1)
        homozygot2_observed_method2 = check_is_number(list_output_cur_homozygot_method2, 1, 1)
        homozygot_alt_method2 = check_is_number(list_output_cur_homozygot_method2, 2, 1) + check_is_number(
            list_output_cur_homozygot_method2, 3, 1) + check_is_number(
            list_output_cur_homozygot_method2, 4, 1)
        heterozygot_alt_method2 = check_is_number(list_output_cur_heterozygot_method2, 1, 1) + check_is_number(
            list_output_cur_heterozygot_method2, 2, 1)

        # Insert the haplotype frequency values into the return DataFrame
        df_to_analyze_return.iloc[len(df_to_analyze) + 9, idxCol] = homozygot1_observed_method2
        df_to_analyze_return.iloc[len(df_to_analyze) + 10, idxCol] = heterozygot_observed_method2
        df_to_analyze_return.iloc[len(df_to_analyze) + 11, idxCol] = homozygot2_observed_method2
        df_to_analyze_return.iloc[len(df_to_analyze) + 12, idxCol] = homozygot_alt_method2
        df_to_analyze_return.iloc[len(df_to_analyze) + 13, idxCol] = heterozygot_alt_method2

        curCategory_method2 = ''
        curCategory_method2 += add_if_not_empty(check_is_empty(list_output_cur_homozygot_method2, 0))
        curCategory_method2 += add_if_not_empty(check_is_empty(list_output_cur_heterozygot_method2, 0))
        curCategory_method2 += add_if_not_empty(check_is_empty(list_output_cur_homozygot_method2, 1))
        #
        frequency_method2 = calculate_frequency(list_output_cur_homozygot_method2, list_output_cur_heterozygot_method2)

        # # Remove the leading backslach
        if curCategory_method2 and curCategory_method2[0] == '/':
            curCategory_method2 = curCategory_method2[1:]
        #
        items_method2 = iter(frequency_method2.items())
        first_item_method2 = next(items_method2, None)
        second_item_method2 = next(items_method2, None)
        third_item_method2 = next(items_method2, None)
        first_item_method2 = check_entry(first_item_method2)
        second_item_method2 = check_entry(second_item_method2)
        third_item_method2 = check_entry(third_item_method2)

        first_item_method2 = first_item_method2[1]
        second_item_method2 = second_item_method2[1] + third_item_method2[1]
        # #
        df_to_analyze_return.iloc[len(df_to_analyze) + 14, idxCol] = first_item_method2
        df_to_analyze_return.iloc[len(df_to_analyze) + 15, idxCol] = second_item_method2
        df_to_analyze_return.iloc[len(df_to_analyze) + 16, idxCol] = curCategory_method2
        ################################################################################################################
        # Hardy Weinberg (uncomment if required)
        ################################################################################################################
        # Hardy Weinberg for method 1:
        homozygot1_expected_method1 = total_method1 * (first_item_method1 ** 2)
        heterozygot_expected_method1 = total_method1 * 2 * first_item_method1 * second_item_method1
        homozygot2_expected_method1 = total_method1 * (second_item_method1 ** 2)
        #
        # Hardy Weinberg for method 2:
        homozygot1_expected_method2 = total_method2 * (first_item_method2 ** 2)
        heterozygot_expected_method2 = total_method2 * 2 * first_item_method2 * second_item_method2
        homozygot2_expected_method2 = total_method2 * (second_item_method2 ** 2)
        #
        # # Create lists of observed and expected frequencies
        observed_method1 = [homozygot1_observed_method1, heterozygot_observed_method1 + heterozygot_alt_method1,
                             homozygot2_observed_method1 + homozygot_alt_method1]
        expected_method1 = [homozygot1_expected_method1, heterozygot_expected_method1, homozygot2_expected_method1]
        expected_method1_rounded = rounded_expected_values(expected_method1, sum(observed_method1))
        #
        if hardy_weinberg ==1:
            df_to_analyze_return.iloc[len(df_to_analyze) + 19, idxCol] = expected_method1_rounded[0]
            df_to_analyze_return.iloc[len(df_to_analyze) + 20, idxCol] = expected_method1_rounded[1]
            df_to_analyze_return.iloc[len(df_to_analyze) + 21, idxCol] = expected_method1_rounded[2]

            # Perform Chi-squared test
            chi2_method1, p_method1 = sp.stats.chisquare(observed_method1, f_exp=expected_method1_rounded)

            if np.isnan(p_method1):
                p_method1 = 1
            #
            # # Create lists of observed and expected frequencies
            observed_method2 = [homozygot1_observed_method2, heterozygot_observed_method2 + heterozygot_alt_method2,
                                 homozygot2_observed_method2 + homozygot_alt_method2]
            expected_method2 = [homozygot1_expected_method2, heterozygot_expected_method2, homozygot2_expected_method2]
            expected_method2_rounded = rounded_expected_values(expected_method2, sum(observed_method2))
            #
            df_to_analyze_return.iloc[len(df_to_analyze) + 24, idxCol] = expected_method2_rounded[0]
            df_to_analyze_return.iloc[len(df_to_analyze) + 25, idxCol] = expected_method2_rounded[1]
            df_to_analyze_return.iloc[len(df_to_analyze) + 26, idxCol] = expected_method2_rounded[2]
            #
            # # Perform Chi-squared test
            chi2_method2, p_method2 = sp.stats.chisquare(observed_method2, f_exp=expected_method2_rounded)
            if np.isnan(p_method2):
                p_method2 = 1
            #
            df_to_analyze_return.iloc[len(df_to_analyze) + 28, idxCol] = np.round(p_method1, 4)
            df_to_analyze_return.iloc[len(df_to_analyze) + 29, idxCol] = np.round(p_method2, 4)

    df_to_analyze_return.replace('nan', "", inplace=True)
    return df_to_analyze_return


def merge_dataframes(overall_haplotypes, output_df):
    # Prepare for merge
    overall_haplotypes.index.name = 'Combinations Overall'
    output_df = output_df[output_df.iloc[:, 0].notnull() & (output_df.iloc[:, 0] != '')]
    overall_haplotypes = overall_haplotypes.sort_values(
        by=[overall_haplotypes.columns[-2], overall_haplotypes.columns[-1]], ascending=[False, False])
    output_df = output_df.sort_values(by=[output_df.columns[-1], output_df.columns[-2]], ascending=[False, True])
    overall_haplotypes.reset_index(inplace=True)
    output_df = output_df.reset_index(drop=True)

    # Move the 'Method' column to the end of output_df
    cols = list(output_df.columns)
    cols.append(cols.pop(cols.index('Method')))
    output_df = output_df[cols]

    # Remove the 'Method' column from overall_haplotypes
    overall_haplotypes = overall_haplotypes.drop(columns=['Method'])

    # Merge overall haplotypes and output dataframes
    output_df_merged = pd.concat([overall_haplotypes, output_df], axis=1)
    output_df_merged.set_index('Combinations Overall', inplace=True)

    return output_df_merged


def create_row_summary(cur_row, positions, cur_amplicon):
    if '_' in cur_row.name and not 'QUAL' in cur_row.name:
        # Initialize an empty list for results
        results = []
        for position, value in zip(positions, cur_row):
            if not pd.isna(value) and 'row summary' not in str(position) and value != '':
                results.append(f"{int(position)}{value}")

        # Join the results with a space
        result = ' '.join(results)

        # Return '' if the result is empty, else return the result
        return '' if not result else result


def row_summary(df_to_analyze, cur_amplicon):
    counter = -1
    ref_positions = list(map(int, df_to_analyze.columns))
    cds_positions = list(map(int, df_to_analyze.loc['variant position by CDS numeration']))
    ncbi_positions = list(map(int, df_to_analyze.loc['variant position by NCBI numeration']))

    for idx, row in df_to_analyze.iterrows():
        counter = counter + 1
        # Skip the rows that provide position information
        if 'numeration' in idx:
            continue

        if '_H1_' in idx:
            # Find the corresponding H2 row
            h2_row = df_to_analyze.loc[idx.replace('_H1_', '_H2_')]

            for pos_set, col_name in zip([ref_positions, cds_positions, ncbi_positions], ['ref.', 'CDS', 'NCBI']):
                # Determine the order based on sequencing direction
                order = pos_set[::-1] if cur_amplicon['Sequencing Direction'] == 'backward' else pos_set

                haplotype_values = []
                for pos, val1 in zip(order, row):
                    val2 = h2_row[df_to_analyze.columns[order.index(pos)]]
                    if not pd.isna(val1):
                        if '/' in str(val1):
                            haplotype_values.append(f"{pos}{val1}")
                        else:
                            haplotype_values.append(f"{pos}{val1}|{val2}")

                haplotype_result = ' '.join(haplotype_values)
                df_to_analyze.at[idx, f'haplotype - {col_name} seq. numeration'] = haplotype_result

        if '_H' in idx:
            for pos_set, col_name in zip([ref_positions, cds_positions, ncbi_positions], ['ref.', 'CDS', 'NCBI']):
                # Determine the order based on sequencing direction
                order = pos_set[::-1] if cur_amplicon['Sequencing Direction'] == 'backward' else pos_set
                cur_row_summary = create_row_summary(row, order, cur_amplicon)
                df_to_analyze.at[idx, f'row summary - {col_name} seq. numeration'] = cur_row_summary

    new_columns_order = [
        'row summary - ref. seq. numeration',
        'haplotype - ref. seq. numeration',
        'row summary - CDS seq. numeration',
        'haplotype - CDS seq. numeration',
        'row summary - NCBI seq. numeration',
        'haplotype - NCBI seq. numeration'
    ]

    # Add the remaining columns to the list
    remaining_columns = [col for col in df_to_analyze.columns if col not in new_columns_order]
    all_columns_ordered = remaining_columns + new_columns_order

    # Reorder the DataFrame columns
    df_to_analyze = df_to_analyze[all_columns_ordered]

    return df_to_analyze


def sort_entry(entry):
    # Check if there is a backslash in the entry
    if '/' in entry:
        # Split the entry at the backslash
        left_side, right_side = entry.split('/')
        # Sort both sides and join them back with a backslash
        sorted_entry = '/'.join([''.join(sorted(left_side)), ''.join(sorted(right_side))])
        return sorted_entry
    # If there is no backslash, return the entry as is
    return entry


def analyze_occurrence(df_to_analyze, methods=None):
    # Define a helper function to extract number from the string
    def extract_number(s):

        if 'x' in str(s):
            return int(s.split('x')[0].strip())
        return None

    # Function to shift non-empty cells to the top
    def shift_cells(column):
        # Filter out empty strings and reset the index without adding the old index as a column
        return column[column != ''].reset_index(drop=True)

    # Sorting function for backslash-separated values in df_H1
    def sort_backslash_separated(entry):
        if '/' in entry:
            return '/'.join(sorted(entry.split('/')))
        return entry

    def process_dataframe(df):
        df_H1 = df[df.index.str.contains('H1')].copy()
        df_H2 = df[df.index.str.contains('H2')].copy()

        df_H1 = df_H1.fillna('')
        df_H2 = df_H2.fillna('')

        df_H1.index = df_H1.index.str.replace('_H1_', '_')
        df_H2.index = df_H2.index.str.replace('_H2_', '_')

        df_H1 = df_H1.sort_index()
        df_H2 = df_H2.sort_index()

        assert (df_H1.index == df_H2.index).all()

        df_H1 = df_H1.applymap(sort_backslash_separated)

        df_joined = df_H1.where(df_H1.applymap(lambda x: '/' in x), df_H1 + df_H2)

        exclude_columns = [col for col in df_joined.columns if "numeration" in str(col)]
        include_columns = [col for col in df_joined.columns if col not in exclude_columns]

        # remove if greater than 6 chars long
        for col in include_columns:
            # Remove backslash from cells that have less than 6 characters
            df_joined[col] = df_joined[col].apply(
                lambda x: x.replace('/', '') if isinstance(x, str) and len(x) < 6 else x)

            # Set to empty string if cell length is more than 6 characters
            df_joined.loc[df_joined[col].astype(str).str.len() > 6, col] = ''

        df_joined_included = df_joined[include_columns].applymap(sort_entry)

        # Count occurrences
        df_counts = df_joined_included.apply(lambda x: x.value_counts())
        df_counts['sum'] = df_counts.sum(axis=1)
        df_counts.sort_values(by='sum', ascending=False, inplace=True)
        df_counts.drop('sum', axis=1, inplace=True)

        df_counts_reset_included = df_counts.reset_index()

        df_combined = df_counts_reset_included.applymap(lambda x: '' if pd.isna(x) else str(x))

        for i in df_combined.index:
            for j in df_combined.columns:
                if j != 'index':
                    if df_combined.at[i, j] != '' and df_combined.at[i, 'index'] != '':
                        df_combined.at[i, j] = str(int(float(df_combined.at[i, j]))) + 'x ' + df_combined.at[i, 'index']
                    else:
                        df_combined.at[i, j] = ''

        df_combined = df_combined.drop(columns=['index'])
        df_combined.columns = df_combined.columns.astype(str)
        # Apply sorting logic to each column
        for col in df_combined.columns:
            # Extract numbers and get the sorted index
            df_combined = df_combined.apply(shift_cells)
            sorted_index = df_combined[col].map(extract_number).sort_values(ascending=False).index
            df_combined[col] = df_combined[col].reindex(sorted_index).reset_index(drop=True).fillna("")
        df_combined = df_combined[df_combined.apply(lambda row: row.str.strip().any(), axis=1)]

        return df_combined

    df_total = process_dataframe(df_to_analyze)
    # Create a dictionary mapping the current index to the new names
    new_index = {i: f"Total allele frequency - allele {i + 1}" for i in range(len(df_total))}

    # Rename the index using the new_index dictionary
    df_total.rename(index=new_index, inplace=True)

    results = [df_to_analyze, df_total]

    if methods:
        method1, method2 = methods
        df_method1 = df_to_analyze[df_to_analyze.index.str.contains(method1)].copy()
        df_method2 = df_to_analyze[df_to_analyze.index.str.contains(method2)].copy()

        df_method1_result = process_dataframe(df_method1)
        df_method2_result = process_dataframe(df_method2)

        # Create a dictionary mapping the current index to the new names
        new_index_method1 = {i: method1 + f" total allele frequency - allele {i + 1}" for i in range(len(df_method1_result))}
        df_method1_result.rename(index=new_index_method1, inplace=True)

        new_index_method2 = {i: method2 + f" total allele frequency - allele {i + 1}" for i in range(len(df_method2_result))}
        df_method2_result.rename(index=new_index_method2, inplace=True)

        results.extend([df_method1_result, df_method2_result])

    results[0].columns = results[0].columns.astype(str)
    df_return = pd.concat(results, axis=0)

    return df_return

def extract_combinations_alleles(df_to_analyse, methods, cur_amplicon):
    # Function to sort values with '/' by alphabet
    def sort_value(value):
        if pd.notna(value) and len(value) == 3 and '/' in value:
            chars = value.split('/')
            return '/'.join(sorted(chars))
        return value

    # Function to combine values from H1 and H2
    def combine_values(value1, value2):
        if pd.isna(value1) or pd.isna(value2):
            return np.nan
        value1, value2 = str(value1), str(value2)
        if value2 < value1:
            value1, value2 = value2, value1
        if value1 == value2 and value1 != '':
            if '/' in value1:
                return value1.replace('/', '/')
            else:
                return value1 + '|' + value2
        else:
            return value1 + '|' + value2


    # Transpose the DataFrame for easier processing
    df_to_analyse = df_to_analyse.T

    # Create a new DataFrame for the results
    results = pd.DataFrame(columns=["alleles by ref.seq. numeration", "alleles by CDS numeration", "alleles by NCBI numeration",
                                    f"Counter {methods[0]}", f"Counter {methods[1]}"])

    # Get the list of positions
    positions = df_to_analyse.index.tolist()

    # Reverse the positions if the sequencing direction is 'backward'
    if cur_amplicon['Sequencing Direction'] == 'backward':
        positions = positions[::-1]

    # Iterate through the columns of the transposed DataFrame
    for col in df_to_analyse.columns:
        combinations_seen = set()  # Set to track unique combinations

        # Iterate through each method
        for method in methods:
            if f"_H1_{method}" in col:
                base_col_name = col.replace(f"_H1_{method}", "")
                h2_col = f"{base_col_name}_H2_{method}"

                # Combine H1 and H2 values
                cur_combination = df_to_analyse.apply(lambda row: combine_values(row[col], row[h2_col]), axis=1)
                cur_combination = cur_combination.dropna().apply(sort_value)

                # Determine the numerations
                items = list(cur_combination.items())
                if cur_amplicon['Sequencing Direction'] == 'backward':
                    haplotype_ref_numeration = " ".join([f"{int(index)}{value}" for index, value in reversed(items)])
                    haplotype_cds_numeration = " ".join([f"{df_to_analyse['variant position by CDS numeration'].loc[index]}{value}" for index, value in reversed(items) if index in df_to_analyse.index])
                    haplotype_ncbi_numeration = " ".join([f"{df_to_analyse['variant position by NCBI numeration'].loc[index]}{value}" for index, value in reversed(items) if index in df_to_analyse.index])
                else:
                    haplotype_ref_numeration = " ".join([f"{int(index)}{value}" for index, value in items])
                    haplotype_cds_numeration = " ".join([f"{df_to_analyse['variant position by CDS numeration'].loc[index]}{value}" for index, value in items if index in df_to_analyse.index])
                    haplotype_ncbi_numeration = " ".join([f"{df_to_analyse['variant position by NCBI numeration'].loc[index]}{value}" for index, value in items if index in df_to_analyse.index])

                if haplotype_ref_numeration != '':
                    # Check if the combination is new or seen before
                    if haplotype_ref_numeration not in results["alleles by ref.seq. numeration"].values:
                        new_row = pd.DataFrame({"alleles by ref.seq. numeration": [haplotype_ref_numeration],
                                                "alleles by CDS numeration": [haplotype_cds_numeration],
                                                "alleles by NCBI numeration": [haplotype_ncbi_numeration],
                                                f"Counter {methods[0]}": [0],
                                                f"Counter {methods[1]}": [0]})
                        results = pd.concat([results, new_row], ignore_index=True)

                    # Update the counters
                    results.loc[(results["alleles by ref.seq. numeration"] == haplotype_ref_numeration) &
                                (results["alleles by CDS numeration"] == haplotype_cds_numeration) &
                                (results["alleles by NCBI numeration"] == haplotype_ncbi_numeration),
                                f"Counter {method}"] += 1

    # Sort results by the first method's counter
    results = results.sort_values(by=f"Counter {methods[0]}", ascending=False).reset_index(drop=True)
    return results

def extract_combinations_haplotypes(df_to_analyse, methods, cur_amplicon):
    # Function to sort values with '/' by alphabet
    def sort_value(value):
        if pd.notna(value) and len(value) == 3 and '/' in value:
            chars = value.split('/')
            return '/'.join(sorted(chars))
        return value

    # Transpose the DataFrame for easier processing
    df_to_analyse = df_to_analyse.T

    # Create a new DataFrame for the results
    results = pd.DataFrame(columns=["haplotypes by ref.seq. numeration", "haplotypes by CDS numeration", "haplotypes by NCBI numeration",
                                    f"Counter {methods[0]}", f"Counter {methods[1]}"])

    # Get the list of positions
    positions = df_to_analyse.index.tolist()

    # Reverse the positions if the sequencing direction is 'backward'
    if cur_amplicon['Sequencing Direction'] == 'backward':
        positions = positions[::-1]

    # Iterate through the columns of the transposed DataFrame
    for col in df_to_analyse.columns:
        combinations_seen = set()  # Set to track unique combinations

        # Iterate through each method
        for method in methods:
            if f"_H1_{method}" in col:
                base_col_name = col.replace(f"_H1_{method}", "")
                h2_col = f"{base_col_name}_H2_{method}"

                # Combine H1 and H2 values
                cur_combination_h1 = df_to_analyse.apply(lambda row: row[col], axis=1)
                cur_combination_h2 = df_to_analyse.apply(lambda row: row[h2_col], axis=1)

                # Determine the numerations
                items_h1 = list(cur_combination_h1.dropna().items())
                items_h2 = list(cur_combination_h2.dropna().items())

                if cur_amplicon['Sequencing Direction'] == 'backward':
                    haplotype_ref_numeration_h1 = " ".join([f"{int(index)}{value}" for index, value in reversed(items_h1)])
                    haplotype_ref_numeration_h2 = " ".join([f"{int(index)}{value}" for index, value in reversed(items_h2)])
                    haplotype_cds_numeration_h1 = " ".join([f"{df_to_analyse['variant position by CDS numeration'].loc[index]}{value}" for index, value in reversed(items_h1) if index in df_to_analyse.index])
                    haplotype_cds_numeration_h2 = " ".join([f"{df_to_analyse['variant position by CDS numeration'].loc[index]}{value}" for index, value in reversed(items_h2) if index in df_to_analyse.index])
                    haplotype_ncbi_numeration_h1 = " ".join([f"{df_to_analyse['variant position by NCBI numeration'].loc[index]}{value}" for index, value in reversed(items_h1) if index in df_to_analyse.index])
                    haplotype_ncbi_numeration_h2 = " ".join([f"{df_to_analyse['variant position by NCBI numeration'].loc[index]}{value}" for index, value in reversed(items_h2) if index in df_to_analyse.index])
                else:
                    haplotype_ref_numeration_h1 = " ".join([f"{int(index)}{value}" for index, value in items_h1])
                    haplotype_ref_numeration_h2 = " ".join([f"{int(index)}{value}" for index, value in items_h2])
                    haplotype_cds_numeration_h1 = " ".join([f"{df_to_analyse['variant position by CDS numeration'].loc[index]}{value}" for index, value in items_h1 if index in df_to_analyse.index])
                    haplotype_cds_numeration_h2 = " ".join([f"{df_to_analyse['variant position by CDS numeration'].loc[index]}{value}" for index, value in items_h2 if index in df_to_analyse.index])
                    haplotype_ncbi_numeration_h1 = " ".join([f"{df_to_analyse['variant position by NCBI numeration'].loc[index]}{value}" for index, value in items_h1 if index in df_to_analyse.index])
                    haplotype_ncbi_numeration_h2 = " ".join([f"{df_to_analyse['variant position by NCBI numeration'].loc[index]}{value}" for index, value in items_h2 if index in df_to_analyse.index])

                # Combine results for H1 and H2
                haplotype_ref_numeration = haplotype_ref_numeration_h1 + " " + haplotype_ref_numeration_h2
                haplotype_cds_numeration = haplotype_cds_numeration_h1 + " " + haplotype_cds_numeration_h2
                haplotype_ncbi_numeration = haplotype_ncbi_numeration_h1 + " " + haplotype_ncbi_numeration_h2

                if haplotype_ref_numeration != '':
                    # Check if the combination is new or seen before
                    if haplotype_ref_numeration not in results["haplotypes by ref.seq. numeration"].values:
                        new_row = pd.DataFrame({"haplotypes by ref.seq. numeration": [haplotype_ref_numeration],
                                                "haplotypes by CDS numeration": [haplotype_cds_numeration],
                                                "haplotypes by NCBI numeration": [haplotype_ncbi_numeration],
                                                f"Counter {methods[0]}": [0],
                                                f"Counter {methods[1]}": [0]})
                        results = pd.concat([results, new_row], ignore_index=True)

                    # Update the counters
                    counter_increment = 2 if cur_combination_h1.equals(cur_combination_h2) else 1
                    results.loc[(results["haplotypes by ref.seq. numeration"] == haplotype_ref_numeration) &
                                (results["haplotypes by CDS numeration"] == haplotype_cds_numeration) &
                                (results["haplotypes by NCBI numeration"] == haplotype_ncbi_numeration),
                                f"Counter {method}"] += counter_increment

    results = results[results["haplotypes by ref.seq. numeration"].str.strip() != '']

    # Sort results by the first method's counter
    results = results.sort_values(by=f"Counter {methods[0]}", ascending=False).reset_index(drop=True)
    return results

# Inserts a row into the given DataFrame at the first position.
def insert_row(df_to_analyze, new_row_values, row_label):
    # Create a copy of the input DataFrame to avoid modifying the original
    df_to_return = df_to_analyze.copy()

    # Create a new DataFrame from new_row_values with the same columns as df_to_analyze
    # Use the row_label as the index for this new row
    new_row = pd.DataFrame([new_row_values], columns=df_to_analyze.columns, index=[row_label])

    # Concatenate the new row DataFrame with the copied DataFrame
    # By default, pd.concat appends the new row DataFrame to the top
    df_to_return = pd.concat([new_row, df_to_return])

    # Return the DataFrame with the inserted row
    return df_to_return


# Insert the reference sequence into the given DataFrame, the reference sequence is added as the first row
def insert_reference_sequence(df_to_analyze, reference_sequence):
    reference_row = [reference_sequence.get(int(float(col)), np.nan) if not isinstance(col, str) else col for col in
                     df_to_analyze.columns]
    return insert_row(df_to_analyze, reference_row, 'reference sequence')


# Process reference sequence
def process_sequence(file_path, amplicon_data):
    # Open and read the file
    with open(file_path, 'r') as f:
        sequences = f.read()

    # Merge all sequences into one string
    merged_sequence = ''.join(sequences)

    # Remove newline characters
    merged_sequence = merged_sequence.replace('\n', '')

    # Prepare an empty dictionary to hold sequence information
    sequence_dict = {}

    # Extract the keys (amplicon names) from amplicon_data dictionary
    amplicons = list(amplicon_data.keys())

    # Iteratively find each amplicon's sequence
    for i in range(len(amplicons) - 1):
        # Compile a regular expression pattern
        pattern = re.compile(amplicons[i] + '(.*?)' + amplicons[i + 1])

        # Search for the pattern in the merged_sequence
        match = pattern.search(merged_sequence)

        # If a match is found
        if match:
            # Extract the sequence
            sequence = match.group(1)

            # Filter out characters other than 'ACTG'
            sequence = ''.join([ch for ch in sequence if ch in 'ACTG'])

            # Create a dictionary with position as key and corresponding character as value
            position_dict = {pos + 1: char for pos, char in enumerate(sequence)}

            # Add this dictionary to sequence_dict with amplicon name as the key
            sequence_dict[amplicons[i]] = position_dict

    # Handle the last amplicon separately
    pattern = re.compile(amplicons[-1] + '(.*?)$')
    match = pattern.search(merged_sequence)

    if match:
        sequence = match.group(1)
        sequence = ''.join([ch for ch in sequence if ch in 'ACTG'])
        position_dict = {pos + 1: char for pos, char in enumerate(sequence)}
        sequence_dict[amplicons[-1]] = position_dict

    # Return the final sequence_dict containing information for all amplicons
    return sequence_dict


# Insert the NCBI (National Center for Biotechnology Information) information to the DataFrame based on the amplicon value.
def insert_ncbi(df_to_analyze, cur_amplicon_data):
    offset = cur_amplicon_data['Offset NCBI']
    reverse = cur_amplicon_data['Sequencing Direction'] == 'backward'
    ncbi_row = [offset - name if not isinstance(name, str) else name for name in df_to_analyze.columns]

    # If sequencing direction is 'backward', reverse the new row except for the 'NCBI' column
    if reverse:
        ncbi_row = df_to_analyze.columns[::-1]

    # Use the existing insert_row function to insert the new row
    df_to_analyze = insert_row(df_to_analyze, ncbi_row, 'variant position by NCBI numeration')

    return df_to_analyze


# Insert the CDS (Coding DNA Sequence) information to the DataFrame based on the amplicon value.
def insert_cds(df_to_analyze, cur_amplicon_data):
    offset = cur_amplicon_data['Offset CDS']
    cds_row = [name - offset if isinstance(name, int) else name for name in df_to_analyze.iloc[0]]
    return insert_row(df_to_analyze, cds_row, 'variant position by CDS numeration')


# Insert empty rows at the end of the DataFrame and assign specific row names.
def insert_empty_rows_method1_method2(df_to_analyze, methods, num_rows, hardy_weinberg):
    method1 = methods[0]
    method2 = methods[1]

    # Append NaN rows to the DataFrame
    for _ in range(num_rows):
        df_to_analyze.loc[len(df_to_analyze.index)] = np.nan

    if hardy_weinberg == 1:
        row_names = {
            idx: name for idx, name in zip(
                df_to_analyze.index[-(num_rows):],
                [
                    '', method1 + ' observed homozygous 1', method1 + ' observed heterozygous',
                    method1 + ' observed homozygous 2', method1 + ' observed homozygous alt',
                    method1 + ' observed heterozygous alt', method1 + ' frequency 1', method1 + ' frequency 2',
                    method1 + ' category', method2 + ' observed homozygous 1', method2 + ' observed heterozygous',
                    method2 + ' observed homozygous 2', method2 + ' observed homozygous alt',
                    method2 + ' observed heterozygous alt', method2 + ' frequency 1', method2 + ' frequency 2',
                    method2 + ' category', '',
                    method1 + ' expected', method1 + ' expected homozygous 1',
                    method1 + ' expected heterozygous', method1 + ' expected homozygous 2', '',
                    method2 + ' expected', method2 + ' expected homozygous 1', method2 + ' expected heterozygous',
                    method2 + ' expected homozygous 2', '', method1 + ' Hardy Weinberg', method2 + ' Hardy Weinberg', ''
                ]
            )
        }
    else:
        row_names = {
            idx: name for idx, name in zip(
                df_to_analyze.index[-(num_rows):],
                [
                    '', method1 + ' observed homozygous 1', method1 + ' observed heterozygous',
                        method1 + ' observed homozygous 2', method1 + ' observed homozygous alt',
                        method1 + ' observed heterozygous alt', method1 + ' frequency 1', method1 + ' frequency 2',
                        method1 + ' category', method2 + ' observed homozygous 1', method2 + ' observed heterozygous',
                        method2 + ' observed homozygous 2', method2 + ' observed homozygous alt',
                        method2 + ' observed heterozygous alt', method2 + ' frequency 1', method2 + ' frequency 2',
                        method2 + ' category', '',
                ]
            )
        }

    # Rename the DataFrame index based on the row_names dictionary
    df_to_analyze = df_to_analyze.rename(index=row_names)
    return df_to_analyze

def annotate_chi_square(cur_result: pd.DataFrame) -> pd.DataFrame:
    # Create a mask for rows containing "Hardy Weinberg" in the row name, safely handling NaN values
    mask = cur_result.index.str.contains("Hardy Weinberg", case=False, na=False)
    chi_square_rows = cur_result[mask]

    # Exclude columns with "row summary" in their name
    columns_to_include = [col for col in cur_result.columns if "row summary" not in col]

    for idx, row in chi_square_rows.iterrows():
        # Filter the row based on the included columns
        row_filtered = row[columns_to_include]

        # Consider empty strings as non-significant
        significant_count = sum((row_filtered < 0.05) & (row_filtered != ''))
        total_values = len(row_filtered)
        not_significant_count = total_values - significant_count

        # Generate the counter string
        counter_string = f"{not_significant_count} of {total_values} are not significant"

        # If "Summary Hardy Weinberg" column doesn't exist, create it right before "row summary - ref. seq. numeration"
        if "Summary Hardy Weinberg" not in cur_result.columns:
            target_idx = cur_result.columns.get_loc("row summary - ref. seq. numeration")
            cur_result.insert(target_idx, "Summary Hardy Weinberg", None)

        # Update the "Summary Hardy Weinberg" column with the counter string
        cur_result.loc[idx, "Summary Hardy Weinberg"] = counter_string

    return cur_result

def sort_excel_sheets(sheet_names):
    # Separate main sheets and haplotypes by methods sheets
    main_sheets = [name for name in sheet_names if not name.endswith('_by_methods')]
    haplotype_sheets = [name for name in sheet_names if name.endswith('_haplotypes_by_methods')]
    allele_sheets = [name for name in sheet_names if name.endswith('_alleles_by_methods')]

    # Sort both lists individually
    main_sheets.sort()
    haplotype_sheets.sort()
    allele_sheets.sort()

    # Combine the sorted lists
    sorted_sheets = main_sheets + haplotype_sheets + allele_sheets

    return sorted_sheets

def analysis_tool(amplicon_data, methods, file_name_merged, path_output_excel, path_reference_sequence,
                  qual_filter_method1, dp_filter_method1, qual_filter_method2, dp_filter_method2):
    # Get the current date
    today = date.today()

    # Define the name of the output excel file
    file_name_analysed = path_output_excel + '/' + str(today) + '_step_3_dataset_analyzed.xlsx'

    # Read data from input Excel file and reference sequence file
    amplicon_names = list(amplicon_data.keys())

    total_result = [pd.read_excel(file_name_merged, sheet_name=sheet) for sheet in amplicon_names]
    ref_sequences = process_sequence(path_reference_sequence, amplicon_data)

    # Prepare dictionary to hold DataFrames for writing later
    output_sheets = {}

    for i, cur_result in enumerate(total_result):
        cur_result.rename(columns={'Unnamed: 0': 'variant position by ref. seq. numeration'}, inplace=True)
        cur_result.set_index('variant position by ref. seq. numeration', inplace=True)
        # Set the title for the row names
        cur_result.index.name = 'variant position by ref. seq. numeration'

        # Get current amplicon
        cur_amplicon_data = amplicon_data[amplicon_names[i]]

        # Get current reference sequence
        cur_ref_sequence = ref_sequences[amplicon_names[i]]

        # Insert NCBI and CDS information and reference sequence
        cur_result = insert_reference_sequence(cur_result, cur_ref_sequence)
        cur_result = insert_ncbi(cur_result, cur_amplicon_data)
        cur_result = insert_cds(cur_result, cur_amplicon_data)

        # Filter data based on QUAL and DP for method1 and method2
        cur_result = filter_data_qual_dp(cur_result, methods, qual_filter_method1, dp_filter_method1,
                                         qual_filter_method2, dp_filter_method2)

        # Apply the cutoff sequence which cuts the df which are under or over the limits
        lower_limit = amplicon_data[amplicon_names[i]]['Lower Limit']
        upper_limit = amplicon_data[amplicon_names[i]]['Upper Limit']
        cur_result = cut_off_sequence(cur_result, lower_limit, upper_limit)

        # Extract allele and haplotype combinations
        overall_result_alleles = extract_combinations_alleles(cur_result, methods, cur_amplicon_data)
        overall_result_haplotypes = extract_combinations_haplotypes(cur_result, methods, cur_amplicon_data)

        # Analyze occurrence and frequencies
        cur_result = row_summary(cur_result, cur_amplicon_data)
        cur_result = analyze_occurrence(cur_result, methods)
        cur_result = extract_haplotype_frequencies(cur_result, methods)
        cur_result.index.name = "variant position by provided reference sequence"

        # Add results to the output_sheets dictionary
        output_sheets[amplicon_names[i]] = cur_result
        output_sheets[f'{amplicon_names[i]}_haplotypes_by_methods'] = overall_result_haplotypes
        output_sheets[f'{amplicon_names[i]}_alleles_by_methods'] = overall_result_alleles

    # Sort the sheet names
    sorted_sheet_names = sort_excel_sheets(list(output_sheets.keys()))

    with pd.ExcelWriter(file_name_analysed, engine='xlsxwriter') as writer:
        # Write each sheet to the Excel file in the sorted order
        for sheet_name in sorted_sheet_names:
            if sheet_name.endswith('_by_methods'):
                output_sheets[sheet_name].to_excel(writer, sheet_name=sheet_name, index=False)
            else:
                output_sheets[sheet_name].to_excel(writer, sheet_name=sheet_name)

    return file_name_analysed