import numpy as np
import pandas as pd

# Sperm information can be stored in two main formats. I call once of them the
# 'Barcode format', and the other is called 'Tuple format'.
#
# Barcode format is a pandas DataFrame in which each barcode gets its own row.
# This is the format output from the data pipeline.
# It may look like:
# 
# barcode         | fish  | sample | replicate | reads
# ----------------+----------------+-----------+-------
# NONE_NONE_NONE  | fish5 | 1      | 1         | 9595
# NONE_1D+37_NONE | fish5 | 1      | 1         | 5871
# NONE_NONE_3D+65 | fish5 | 1      | 1         | 4091
#
# 'Tuple' format lists the reads from all barcodes in a single row. This is the
# format output by model runs and input by likelihood calculators.
# Example:
#
# fish  | sample | replicate | obs_result         | n_sperm
# ------+--------+-----------+--------------------+--------
# fish5 | 1      | 1         | (9595, 5871, 4091) | 19557

sample_map_int = {1:5, 2:6, 3:7, 4:8, 5:9, 6:10, 7:11, 8:12, 9:13, 10:14,
                  11:16, 12:17, 13:18, 14: 20, 15: 22, 16: 24, 17: 25, 18:26}

sample_map = {1:4.86, 2:5.76, 3:6.73, 4:7.86, 5:8.86, 6:9.93, 7:11.0, 8:12.1, 9:13.2, 10:14.23,
                  11:16.1, 12:17.13, 13:18.4, 14:20.53, 15:22.43, 16:23.6, 17: 24.7, 18:25.6}

def get_barcodes(data):
    """
    Inputs a dataset in barcode format and returns a list of the barcodes
    that appear. Barcodes are listed in increasing order of the first
    sample it appears and decreasing order of reads.
    """

    return data.sort_values(
        by = ['sample', 'reads'],
        ascending = [True, False]
    )['barcode'].unique()

#
#
# This function is currently unused because I made an edited version in
# `convert_experimental_data` that works off of frequencies. Consider
# repurposing or removing.
# def convert_to_tuple_format(data, unique_barcodes):
#     """
#     Load the data corresponding to a given `fish`
#     and format it in a dataframe suitable for use
#     in my model.

#     The returned dataframe has columns `time`, `replicate`,
#     and `obs_result`.
#     """
    
#     barcode_index = {s: i for i, s in enumerate(unique_barcodes)}

#     def create_read_arr(barcode, read_count):
#         """
#         Create an array of length `len(barcode_index)` that contains
#         `read_count` at position `barcode_index[barcode]`.
#         """
#         reads = np.zeros(len(unique_barcodes), dtype = int)
#         reads[barcode_index[barcode]] = read_count
#         return reads

#     data['obs_result'] = data.apply(lambda row: create_read_arr(row['barcode'], row['reads']), axis = 1)

#     # Group by sample and replicate, and sum over all the 'obs_result' arrays.
#     grouped_df = data.groupby(['sample', 'replicate'], as_index=False)['obs_result'].sum()

#     grouped_df['obs_result'] = grouped_df['obs_result'].apply(tuple)
#     grouped_df = grouped_df.rename(columns={'sample': 'time'})

#     # Create a new series, `isempty`, that is `True` if the clone of interest has size 0.
#     # We use this to remove all rows of the dataframe before the clone actually appears.
#     isempty = grouped_df['obs_result'].apply(lambda x: x[0] == 0)
#     smallest_time_with_clone = grouped_df[isempty == False]['time'].min()
#     grouped_df = grouped_df[grouped_df['time'] >= smallest_time_with_clone]

#     grouped_df = grouped_df.reset_index(drop=True)
#     grouped_df['n_sperm'] = grouped_df['obs_result'].apply(sum)
    
#     return grouped_df

def convert_to_tuple_format(data, n_sperm = None):
    """
    Load the data corresponding to a given `fish`
    and format it in a dataframe suitable for use
    in my model.

    The returned dataframe has columns `time`, `replicate`,
    and `obs_result`.
    """
    data = data.copy()

    unique_barcodes = get_barcodes(data)
    barcode_index = {s: i for i, s in enumerate(unique_barcodes)}

    def create_read_arr(barcode, proportion):
        """
        Create an array of length `len(barcode_index)` that contains
        `read_count` at position `barcode_index[barcode]`.
        """

        if n_sperm is None:
            reads = np.zeros(len(barcode_index), dtype = float)
            reads[barcode_index[barcode]] = proportion
        else:
            reads = np.zeros(len(barcode_index), dtype = int)
            reads[barcode_index[barcode]] = round(proportion * n_sperm)
        return reads

    data['obs_result'] = data.apply(lambda row: create_read_arr(row['barcode'], row['proportion']), axis = 1)

    # Group by sample and replicate, and sum over all the 'obs_result' arrays.
    grouped_df = data.groupby(['sample', 'replicate'], as_index=False)['obs_result'].sum()

    grouped_df['obs_result'] = grouped_df['obs_result'].apply(tuple)
    grouped_df = grouped_df.rename(columns={'sample': 'time'})

    # Create a new series, `isempty`, that is `True` if the clone of interest has size 0.
    # We use this to remove all rows of the dataframe before the clone actually appears.
    isempty = grouped_df['obs_result'].apply(lambda x: x[0] == 0)
    smallest_time_with_clone = grouped_df[isempty == False]['time'].min()
    grouped_df = grouped_df[grouped_df['time'] >= smallest_time_with_clone]

    grouped_df = grouped_df.reset_index(drop=True)
    
    if n_sperm is not None:
        grouped_df['n_sperm'] = n_sperm #grouped_df['obs_result'].apply(sum)
    
    return grouped_df

def add_tuples_elementwise(tuples):
    """
    This function inputs a series of tuples, and outputs their elementwise
    sum. All tuples must have the same length.
    
    Example:
    add_tuples_elementwise((1, 2), (3, 4), (5, 6)) --> (9, 12)
    """
    length = len(tuples.iloc[0])
    if not all(len(t) == length for t in tuples):
        raise ValueError("All input tuples must have the same length.")
    
    return tuple(sum(values) for values in zip(*tuples))

# def normalize_tuple(tup, target_sum):
#     """
#     This function inputs a tuple `tup` and divides each entry a constant
#     so that the sum of the entries equals `target_sum`.

#     The purpose of this function is to convert a tuple of sperm 'reads'
#     from each clone into plausible sperm cell counts.
    
#     tup: input tuple
#     count: target total sum
#     """
#     return tuple([round(idx/sum(tup)*target_sum) for idx in tup])

def normalize_tuple(tup, target_sum):
    """
    Scale the values in `tup` so that their sum equals `target_sum`,
    while keeping integer outputs. Uses the largest remainder method
    to fairly distribute rounding errors.
    """

    # Step 1: Compute the scaling factor
    total = sum(tup)
    scaled_values = [(value / total) * target_sum for value in tup]

    # Step 2: Take the floor of each scaled value (convert to int)
    floored_values = [int(val) for val in scaled_values]

    # Step 3: Calculate how much we are short of the target sum
    current_sum = sum(floored_values)
    remainder = target_sum - current_sum

    # Step 4: Compute fractional parts of each scaled value
    fractional_parts = [val - int(val) for val in scaled_values]

    # Step 5: Identify indices with the largest fractional parts
    # These are the ones that should get the extra counts
    indices_sorted_by_fraction = sorted(
        range(len(tup)),
        key=lambda idx: fractional_parts[idx],
        reverse=True
    )

    # Step 6: Distribute the remainder across the top indices
    for idx in indices_sorted_by_fraction[:remainder]:
        floored_values[idx] += 1

    # Step 7: Return the adjusted tuple
    return tuple(floored_values)

def normalize_sperm_counts(data, target_sum):
    """
    This function takes a dataset in Tuple format and normalizes the sperm counts
    to equal a target sum.
    """
    data = data.copy()
    data['obs_result'] = data['obs_result'].apply(normalize_tuple, target_sum=target_sum)
    data['n_sperm'] = data['obs_result'].apply(sum)

    return data

def reorder_tuple(tup, order):
    """Reorders a tuple according to the given order."""
    return tuple(tup[i] for i in order)

def sort_results(data, argsort = None):
    """
    Sort the obs_result (and hid_result, if applicable) column of a Tuple
    formatted dataset.

    If `argsort` is not specified, then the ordering is done so that clones
    are ordered in increasing order of the time they appear, and decreasing
    in the number of reads.

    If `argsort` is specified, then clones are reordered according to the
    indices of `argsort`.
    """
    obs_result = data['obs_result']

    if argsort is None:
        n_clones = len(obs_result.iloc[0])
        samples = []
        reads = []
        
        for i in range(n_clones):
            for l, tup in enumerate(obs_result):
                if tup[i] != 0:
                    samples.append(l)
                    reads.append(tup[i])
                    break
            else:
                # If all entries are zero for this index i
                samples.append(None)
                reads.append(0)
        
        # Create the new dataframe
        df = pd.DataFrame({'sample': samples, 'reads': reads})
        argsort = df.sort_values(
                by = ['sample', 'reads'],
                ascending = [True, False]
            ).index.to_numpy()

    data['obs_result'] = data['obs_result'].apply(lambda tup: reorder_tuple(tup, argsort))
    if 'hid_result' in data.columns:
        data['hid_result'] = data['hid_result'].apply(lambda tup: reorder_tuple(tup, argsort))
    if 'full_result' in data.columns:
        data['full_result'] = data['full_result'].apply(lambda tup: reorder_tuple(tup, argsort))
        
    return argsort

def collapse(data, idx):
    """
    Reduce all but one stem cell clone into a single group
    for direct comparison with a clone of interest.

    Input `results` is a data frame of model results, and `idx`
    is the index of the clone we want to isolate.
    """
    data = data.copy()

    # print(data)

    def combine(t):
        """
        Input a tuple `t` of any length, and return
        a 2-tuple with the first element equal to `t[idx]`
        and second element equal to the sum of the remaining elements.
        """
        selected = t[idx]
        remaining_sum = sum(t[:idx] + t[idx+1:])
        return (selected, remaining_sum)

    data['obs_result'] = data['obs_result'].apply(lambda x: combine(x))
    if 'hid_result' in data.columns:
        data['hid_result'] = data['hid_result'].apply(lambda x: combine(x))
    if 'full_result' in data.columns:
        data['full_result'] = data['full_result'].apply(lambda x: combine(x))
        
    data['clone_freq'] = data['obs_result'].apply(lambda result: result[0]) / data['n_sperm']
        
    return data

def where_unedited(barcode_list):
    """
    Given a list of barcodes (strings), returns the index where the unedited barcode occurs.
    If the unedited barcode does not appear, return `None` instead.
    """
    unedited_barcode = 'NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE'
    where = np.where(barcode_list == unedited_barcode)[0]
    return where[0] if len(where) == 1 else None

def combine_replicates(data):
    """
    Input data in table format and combine the replicates together.

    This function is currently unused.
    """
    agg = (
        data.groupby("time", as_index=False)
          .agg({
              "obs_result": tm.data_processing.add_tuples_elementwise,  # elementwise sum of tuples
              "n_sperm": "sum"                                         # sum across replicates
          })
    )
    
    #Recompute clone_freq from aggregated obs_result
    agg["clone_freq"] = agg["obs_result"].apply(
        lambda x: x[0] / (x[0] + x[1]) if (x[0] + x[1]) > 0 else float("nan")
    )

    return agg

#
# Generating and loading simulation data
#

