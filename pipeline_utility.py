# pipeline_utility.py
# 
# This file contains numerous helper functions for use in the main pipeline
# script, `process_sperm_barcodes.ipynb`. This helps reduce the length of
# the main pipeline and maintain its readability.
# 
# Authors: Connor Shrader & Drew Shapiro

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import networkx as nx
import Levenshtein
import os
from matplotlib.patches import Patch
import re

#
# Reading and writing to files
# 

# def directory(fish, sample, replicate):
#     """
#     This function takes in an identifier for the fish, sample, and replicate,
#     and returns a string representing the location of files for that sample.

#     For example, calling directory('fish1', '2', '3') will return
#     './SpermAnalysis_Muller1/data/viz-Output/fish1_s2_rep3/fish1_s2_rep3'.
#     """
#     filename = f"{fish}_s{sample}_rep{replicate}"
#     return f"./data/{dataset}/{filename}/{filename}"

def tuple_from_file(filepath):
    """
    This function takes in a filepath and returns the fish, sample, and replicate
    that correspond to it. It is, in a sense, the inverse function of directory.

    For example, calling
    directory('./SpermAnalysis_Muller1/data/viz-Output/fish1_s2_rep3/fish1_s2_rep3.allReadCounts')
    will return ('fish1', '2', '3').
    """
    filename = os.path.splitext(os.path.basename(filepath))[0]
    fish, sample, replicate = filename.split('_')
    # fish, sample, replicate = (filepath.split("/")[-1]).split('.')[0].split('_')
    return fish, sample[1:], replicate[3:]

#
# Part C: Barcode clustering
#

def compute_levenshtein_matrix(sequences):
    """
    This function computes the pairwise distances between an array
    of DNA sequences.

    Input: `sequences`: An array of DNA sequences (of length n).
    
    Output: `D`: An (n x n) matrix where D[i, j] is the Levenshtein
                 distance between sequence i and sequence j.
    """
    n = len(sequences)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            d = Levenshtein.distance(sequences[i], sequences[j])
            D[i, j] = D[j, i] = d
    return D

def find_connected_components(D, threshold=1):
    d_inf = D.copy()
    np.fill_diagonal(d_inf, np.inf)
    G = nx.from_numpy_array(d_inf <= threshold)
    return list(nx.connected_components(G))

def summarize_component(i, component_data, distance_matrix):
    """
    This function inputs information about a connected component from
    clustering and outputs a string describing the cluster.

    `i`: The component index.
    `component_data`: All data associated with the cluster.
    `distance_matrix`: The distances between sequences in this cluster.
    """
    output = []
    
    idx = component_data.index
    local_distances = distance_matrix[np.ix_(idx, idx)]
    avg = np.mean(local_distances, axis=0)
    center_idx = np.argmin(avg)

    max_pairwise = int(np.max(np.triu(local_distances)))
    max_to_center = int(np.max(local_distances[center_idx, :]))

    output.append(f"Component {i}\n")
    output.append(f"Max distance between barcodes: {max_pairwise}. Max distance to center: {max_to_center}\n")
    
    for j, row in component_data.iterrows():
        marker = '*' if j == idx[center_idx] else ' '
        output.append(f"{row['reads']:>6}\t{marker}{row['barcode']}.\n")

    return output

def resolve_winner_barcode(component_data, unedited):
    return unedited if unedited in component_data['barcode'].values else component_data.iloc[0]['barcode']

#
# Part E: cleaning contaminated barcodes
#

def identify_questionable(data, unedited):
    """
    This function takes in a DataFrame of reads and identifies which barcodes are questionable.
    It returns a list, called questionable_barcodes, containing the names of barcodes that appear
    in multiple fish.
    """
    barcodes = data[data["barcode"] != unedited].barcode.unique()
    questionable_barcodes = []

    for barcode in barcodes:
        if len(data[data.barcode == barcode].fish.unique()) > 1:
            questionable_barcodes.append(barcode)

    return questionable_barcodes
    
def summarize_barcode(data, barcode):
    """
    This function inputs a DataFrame ('data') of reads, and a barcode. It returns
    a DataFrame ('group_fish') that contains the total number of reads for the
    given barcode in each fish. In addition, this function returns the 'ratio'
    used to identify winner fish.
    """
    single_barcode = data.query('barcode == @barcode').drop("barcode", axis = 1)
    group_fish = single_barcode.groupby(['fish'], observed = True, as_index = False).aggregate({'reads': 'sum'})
    group_fish.sort_values('reads', ascending = False, inplace = True)

    winner_reads = group_fish.reads.iloc[0]
    total_reads = group_fish.reads.sum()

    ratio = winner_reads / (total_reads - winner_reads)

    return group_fish, ratio

def heatmap_barcode(data, barcode):
    """
    This function inputs a DataFrame ('data') of reads and a barcode. It plots a heatmap showing
    the reads of that barcode across all fish, samples, and replicates.
    """
    # data['sample'] = pd.to_numeric(data['sample'])
    group_fish, ratio = summarize_barcode(data, barcode)

    # `idx` is the order that the fish will appear in the heatmap.
    idx = group_fish.index#[::-1]
    
    single_barcode = data[data.barcode == barcode].drop("barcode", axis = 1)
    single_barcode["rep_count"] = 1
    group_samples = single_barcode.groupby(['fish', 'sample'], as_index = False, observed = True) \
                                  .aggregate({'rep_count': 'sum', 'reads': 'sum'})
    
    to_plot = group_samples.pivot(index = 'fish', columns = 'sample', values = 'rep_count') \
                           .fillna(0).iloc[idx]
    
    read_counts = group_samples.pivot(index = "fish", columns = "sample", values = "reads") \
                               .fillna(0).astype(int).iloc[idx]

    # the maximum number of replicates that appear in any fish/sample for this barcode
    max_rep = int(to_plot.values.max())

    # [black, red, green, blue]
    myColors = ((0.0, 0.0, 0.0, 1.0), (0.8, 0.0, 0.0, 1.0), (0.0, 0.8, 0.0, 1.0), (0.0, 0.0, 0.8, 1.0))
    cmap = LinearSegmentedColormap.from_list('Custom', myColors[:max_rep+1], max_rep+1)
    
    fig, ax = plt.subplots()
    ax = sns.heatmap(to_plot, annot = read_counts, annot_kws = {"size": 6}, fmt = 'd', cmap = cmap, vmin = 0, vmax = max_rep)
    # ax.invert_yaxis()
    ax.set_title(f"{barcode}, r={ratio:.2f}", fontsize = 8)
    
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([max_rep/(max_rep+1)*i+max_rep/(2*(max_rep+1)) for i in range(max_rep + 1)])
    colorbar.set_ticklabels([f'{i}' for i in range(max_rep+1)])
    colorbar.set_label("Number of replicates the barcode appeared in")

    return fig, ax

def heatmap_fish(data, fish):
    """
    This function inputs a DataFrame ('data') of reads and a fish. It plots a heatmap showing
    the reads of that barcode across all barcodes, samples, and replicates.
    """
    single_barcode = data.query('fish == @fish').copy().drop("fish", axis = 1)

    # Find appropriate order
    sorted_data = single_barcode.sort_values(by=['sample', 'reads'], ascending=[True, False])
    barcode_order = sorted_data['barcode'].unique()
    
    single_barcode["rep_count"] = 1
    group_samples = single_barcode.groupby(['barcode', 'sample'], as_index = False, observed = True) \
                             .aggregate({'rep_count': 'sum', 'reads': 'sum'})
    
    to_plot = group_samples.pivot(index = 'barcode', columns = 'sample', values = 'rep_count') \
                           .fillna(0).loc[barcode_order]

    # print(single_barcode)
    
    read_counts = group_samples.pivot(index = "barcode", columns = "sample", values = "reads") \
                               .fillna(0).astype(int).loc[barcode_order]

    # the maximum number of replicates that appear in any fish/sample for this barcode
    max_rep = int(to_plot.values.max())

    # [black, red, green, blue]
    myColors = ((0.0, 0.0, 0.0, 1.0), (0.8, 0.0, 0.0, 1.0), (0.0, 0.8, 0.0, 1.0), (0.0, 0.0, 0.8, 1.0))
    cmap = LinearSegmentedColormap.from_list('Custom', myColors[:max_rep+1], max_rep+1)
    
    fig, ax = plt.subplots()
    ax = sns.heatmap(to_plot, annot = read_counts, annot_kws = {"size": 6}, fmt = 'd', cmap = cmap, vmin = 0, vmax = max_rep)
    # ax.invert_yaxis()
    ax.set_title(fish)
    
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([max_rep/(max_rep+1)*i+max_rep/(2*(max_rep+1)) for i in range(max_rep + 1)])
    colorbar.set_ticklabels([f'{i}' for i in range(max_rep+1)])
    colorbar.set_label("Number of replicates the barcode appeared in")

    return fig, ax

def where_unedited(barcode_list):
    """
    Given a list of barcodes (strings), returns the index where the unedited barcode occurs.
    If the unedited barcode does not appear, return `None` instead.
    """
    unedited_barcode = 'NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE'
    where = np.where(barcode_list == unedited_barcode)[0]
    return where[0] if len(where) == 1 else None

# sample_map is a dictionary that maps sample numbers to the actual age
# of the fish during that month.
sample_map = {
    1: 5, 
    2: 6, 
    3: 7, 
    4: 8, 
    5: 9, 
    6: 10, 
    7: 11, 
    8: 12, 
    9: 13, 
    10: 14, 
    11: 16, 
    12: 17, 
    13: 18, 
    14: 20, 
    15: 22, 
    16: 24, 
    17: 25, 
    18: 26
}

def barcodes_to_colors(df):
    sort = df.sort_values(
        by=['sample', 'reads'],
        ascending=[True, False]
    )

    n_barcodes = sort['barcode'].nunique()
    ordered_barcodes = sort['barcode'].unique()
    
    colors = np.array(["#ff4d6f", "#579ea4", "#df7713", "#f9c000", "#86ad34",
              "#5d7298", "#81b28d", "#7e1a2f", "#2d2651", "#c8350d", "#bd777a",
              '#1a5f7a', '#8b4c8f', '#4a9eff', '#ffd93d', '#3dcc8e', '#ff8c42',
              '#e63946', '#a0c334', '#6b4e71', '#00b4d8'])
    unedited_barcode = 'NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE'
    
    repeats = -(-n_barcodes // len(colors))
    colors_extended = np.tile(colors, repeats)[:n_barcodes]

    color_dict = dict()

    i = 0
    for barcode in sort['barcode'].unique():
        if barcode == unedited_barcode:
            color_dict[barcode] = '#aaaaaa'
        else:
            color_dict[barcode] = colors_extended[i]
            i = i + 1

    return ordered_barcodes, color_dict

def plot_muller(data, fish, ax = None, replicate = None, add_legend = False, map_samples = False, color_dict = None, order = None):
    """
    Plot a Muller plot of the barcode data.

    data: The dataset to use (in the pipeline, this will be all_data, pass_data,
        etc.)
    fish: The fish number to be plotted.
    replicate (default None): Which replicate (1 through 3) to plot. If `None`
        is passed, then all three replicates are combined.
    add_legend (defualt False): Whether to add a legend that lists the names
        of barcodes.
    map_samples (default False): If `False`, then the x axis contains the sample
        numbers. If `True`, then the x axis contains the real-life months that
        the samples occurred.
    """
    
    subset = data.query('fish == @fish')
    if replicate is not None:
        subset = subset.query('replicate == @replicate')

    # sorted_data = subset.sort_values(by=['sample', 'reads'], ascending=[True, False])
        
    # barcode_order = sorted_data['barcode'].unique()
    
    # grouped = subset.groupby(['barcode', 'sample'], as_index = False).aggregate({'reads': 'sum'})
    # matrix = grouped.pivot(index = 'barcode', columns = 'sample', values = 'reads')
    # matrix.fillna(0, inplace = True)
    # matrix = matrix.loc[barcode_order]
    # freqs = matrix.div(matrix.sum(axis=0), axis=1)
    
    # Source: https://github.com/BlakeRMills/MoMAColors/blob/main/R/Colors.R
    # colors = np.array(["#ff4d6f", "#579ea4", "#df7713", "#f9c000", "#86ad34",
    #           "#5d7298", "#81b28d", "#7e1a2f", "#2d2651", "#c8350d", "#bd777a",
    #           '#1a5f7a', '#8b4c8f', '#4a9eff', '#ffd93d', '#3dcc8e', '#ff8c42',
    #           '#e63946', '#a0c334', '#6b4e71', '#00b4d8'])
    
    # n_clones = np.shape(freqs)[0]
    # repeats = -(-n_clones // len(colors))
    # colors_extended = np.tile(colors, repeats)[:n_clones]

    if color_dict is None or order is None:
        order, color_dict = barcodes_to_colors(subset)
        

    # barcode_order = sorted_data['barcode'].unique()

    
    grouped = subset.groupby(['barcode', 'sample'], as_index = False).aggregate({'reads': 'sum'})
    matrix = grouped.pivot(index = 'barcode', columns = 'sample', values = 'reads')
    matrix.fillna(0, inplace = True)

    order = [x for x in order if x in matrix.index]
    matrix = matrix.loc[order]
    freqs = matrix.div(matrix.sum(axis=0), axis=1)
    
    colors_extended = [color_dict[bc] for bc in order]
    
    # unedited_index = where_unedited(freqs.index)
    # if unedited_index is not None:
    #     # put in grey in the element corresponding to the unedited array.
    #     colors_extended = np.insert(colors_extended, unedited_index, '#aaaaaa')
    
    #     # remove the last index so that the array maintains the correct length.
    #     colors_extended = colors_extended[:-1]
    
    ts = matrix.columns
    if map_samples:
        ts = [sample_map[t] for t in ts]

    if ax is None:
        fig, ax = plt.subplots()
        
    ax.stackplot(ts, freqs[::-1], colors = colors_extended[::-1]);
    ax.set_xticks(ts);
    ax.set_xlabel("Time")
    ax.set_ylabel("Barcode frequency")
    ax.set_xlim(ts[0], ts[-1])
    ax.set_ylim(0, 1);

    # Create legend handles
    legend_labels = order
    legend_colors = colors_extended
    handles = [Patch(facecolor=c, label=b) for c, b in zip(legend_colors, legend_labels)]
    
    # Add legend to plot
    if add_legend:
        ax.legend(handles=handles[:20], title="Barcode", bbox_to_anchor=(1.05, 1), loc='upper left')

    if ax is None:
        return fig, ax

# Convert barcodes to sequences

def get_sequence(bc: str, ref: str) -> str:
    """
    Get the original sequence from HMID barcode. Applies edits to the reference sequence to obtain a 'clean' version
    of the original sequence that produced the barcode when run through GESTALT.

    :param bc: The HMID barcode
    :param ref: The reference sequence
    :return: The original sequence

    *Example:*\n
    >>>get_sequence('NONE_3D+10_NONE', 'AAACCCGGGTTTAAA')\n
    AAACCCGGGTAA
    """
    delimiters = r'[_&]'
    split = re.split(delimiters, bc)
    split_set = set(split)
    mutation_set = set()

    for edit in split_set:
        if edit == 'NONE':
            continue
        elif 'D' in edit:
            del_length, del_index = map(int, edit.split('D+'))
            mutation_set.add(Mutator(del_index, -del_length, lambda seq, idx, dl=del_length: seq[0:idx] + seq[idx + dl:]))
        elif 'I' in edit:
            ins_length, ins_index, ins_seq = edit.replace('I', '').split('+')
            mutation_set.add(Mutator(int(ins_index), int(ins_length), lambda seq, idx, iseq=ins_seq: seq[0:idx] + iseq + seq[idx:]))
        elif 'S' in edit:
            swp_length, swp_index, swp_seq = edit.replace('S', '').split('+')
            mutation_set.add(Mutator(int(swp_index), 0, lambda seq, idx, sl=int(swp_length), sseq=swp_seq: seq[0:idx] + sseq + seq[idx + sl:]))

    sorted_mutation_set = sorted(mutation_set, key=lambda m: m.start_index)
    index_offset = 0
    for mut in sorted_mutation_set:
        ref = mut.apply(ref, index_offset)
        index_offset += mut.offset_index
    return ref

class Mutator:
    """
    Helper class for applying mutations to a sequence.
    """
    def __init__(self, start_index, offset_index, function):
        self.start_index = start_index
        self.offset_index = offset_index
        self.function = function

    def apply(self, sequence, idx_offset):
        return self.function(sequence, self.start_index + idx_offset)

class BarcodeRemap:
    # The default reference sequence is the zebrafish reference sequence used in the GESTALT pipeline
    zebrafish_reference_sequence = 'AGATCTCGAGCTCaagcttcgGACAGCAGTATCATCGACTATGGagtcGAGAGCGCGCTCGTCGACTATGGagtcGTCAGCAGTACTACTGACGATGGagtcGACAGCAGTGTGTGAGTCTATGGagtcGAGAGCATAGACATCGAGTATGGagtcGACTACAGTCGCTACGACTATGGagtcGACAGAGATATCATGCAGTATGGagtcGACAGCAGTATCTGCTGTCATGGagtcGACTGCACGACAGTCGACTATGGagtcGCGAGCGCTATGAGCGACTATGGGCATGCaattaagtctcagccaccgttaactgaacatgtcaaaacc'

    def __init__(self, reference_sequence=zebrafish_reference_sequence):
        self.reference_sequence = reference_sequence

    def get_all_sequences(self, file_path):
        """
        Finds the 'original' sequences for a set of HMID barcodes derived from the output to the GESTALT pipeline (.allReadCounts or .topReadCounts files).
        Takes the path to a .allReadCounts or .topReadCounts file and returns a pandas dataframe representing the contents of the file, plus a new column
        named 'sequence'. The 'sequence' column contains the original DNA sequences for each barcode.

        :param file_path: The path to the .allReadCounts file
        :return: A pandas dataframe representing the .allReadCounts file with the 'sequence' column added.
        """
        read_count_df = pd.read_csv(file_path, sep='\t')
        read_count_df['sequence'] = read_count_df['event'].apply(lambda x: get_sequence(x, ref=self.reference_sequence))
        return read_count_df