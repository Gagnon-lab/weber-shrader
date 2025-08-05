import numpy as np
import pandas as pd

# Sperm information can be stored in two main formats. I call once of them the
# 'data format', and the other is called 'table format'.
#
# Data format is a pandas DataFrame in which each barcode gets its own row.
# This is the format output from the data pipeline.
# It may look like:
# 
# barcode         | fish  | sample | replicate | reads
# ----------------+----------------+-----------+-------
# NONE_NONE_NONE  | fish5 | 1      | 1         | 9595
# NONE_1D+37_NONE | fish5 | 1      | 1         | 5871
# NONE_NONE_3D+65 | fish5 | 1      | 1         | 4091
#
# Table format lists the reads from each barcode in each row. Example:
#
# fish | sample | replicate | obs_result | n_sperm
# -----+--------+-----------+------------+--------
# fish5 | 1 | 1 | (9595, 5871, 4091) | 19557

def get_barcodes(data):
    """
    Return the barcodes of dataset `data`, in increasing order of the first
    sample it appears and decreasing order of reads.
    """

    return data.sort_values(
        by = ['sample', 'reads'],
        ascending = [True, False]
    )['barcode'].unique()

def convert_data_to_results(data, unique_barcodes):
    """
    Load the data corresponding to a given `fish`
    and format it in a dataframe suitable for use
    in my model.

    The returned dataframe has columns `time`, `replicate`,
    and `obs_result`.
    """
    
    barcode_index = {s: i for i, s in enumerate(unique_barcodes)}

    def create_read_arr(barcode, read_count):
        """
        Create an array of length `len(barcode_index)` that contains
        `read_count` at position `barcode_index[barcode]`.
        """
        reads = np.zeros(len(unique_barcodes), dtype = int)
        reads[barcode_index[barcode]] = read_count
        return reads

    data['obs_result'] = data.apply(lambda row: create_read_arr(row['barcode'], row['reads']), axis = 1)

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
    grouped_df['n_sperm'] = grouped_df['obs_result'].apply(sum)
    
    return grouped_df

def add_tuples_elementwise(*tuples):
    """
    Add multiple tuples elementwise. All tuples must have the same length.
    
    Example:
    add_tuples_elementwise((1, 2), (3, 4), (5, 6)) --> (9, 12)
    """
    length = len(tuples[0])
    if not all(len(t) == length for t in tuples):
        raise ValueError("All input tuples must have the same length.")
    
    return tuple(sum(values) for values in zip(*tuples))

def normalize_tuple(tup, target_sum):
    """
    This function inputs a tuple `tup` and divides each entry a constant
    so that the sum of the entries equals `target_sum`.

    The purpose of this function is to convert a tuple of sperm 'reads'
    from each clone into plausible sperm cell counts.
    
    tup: input tuple
    count: target total sum
    """
    return tuple([round(idx/sum(tup)*target_sum) for idx in tup])

def normalize_sperm_counts(data, target_sum):
    data = data.copy()
    data['obs_result'] = data['obs_result'].apply(normalize_tuple, target_sum=target_sum)
    data['n_sperm'] = data['obs_result'].apply(sum)

    return data

def reorder_tuple(tup, order):
    """Reorders a tuple according to the given order."""
    return tuple(tup[i] for i in order)

def sort_results(data, argsort = None):
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
        
    return data

def where_unedited(barcode_list):
    """
    Given a list of barcodes (strings), returns the index where the unedited barcode occurs.
    If the unedited barcode does not appear, return `None` instead.
    """
    unedited_barcode = 'NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE'
    where = np.where(barcode_list == unedited_barcode)[0]
    return where[0] if len(where) == 1 else None