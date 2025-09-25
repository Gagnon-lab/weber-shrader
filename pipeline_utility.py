# pipeline_utility.py
# 
# This file contains numerous helper functions for use in the main pipeline
# script, `process_sperm_barcodes.ipynb`.
# 
# Author: Connor Shrader

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import networkx as nx
import Levenshtein

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

def summarize_component(i, component_data, distance_matrix, summary):
    idx = component_data.index
    local_distances = distance_matrix[np.ix_(idx, idx)]
    avg = np.mean(local_distances, axis=0)
    center_idx = np.argmin(avg)

    max_pairwise = int(np.max(np.triu(local_distances)))
    max_to_center = int(np.max(local_distances[center_idx, :]))

    summary.write(f"Component {i}: indices {list(idx)}\n")
    summary.write(f"Max distance between barcodes: {max_pairwise}. Max distance to center: {max_to_center}\n")
    
    for j, row in component_data.iterrows():
        marker = '*' if j == idx[center_idx] else ' '
        summary.write(f"{row['reads']:>6}\t{marker}{row['barcode']}.\n")

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
    single_barcode = data[data.barcode == barcode].drop("barcode", axis = 1)
    group_fish = single_barcode.groupby(['fish'], observed = True, as_index = False).aggregate({'reads': 'sum'})
    group_fish.sort_values('reads', ascending = False, inplace = True)

    winner_reads = group_fish.reads.iloc[0]
    total_reads = group_fish.reads.sum()

    ratio = winner_reads / (total_reads - winner_reads)

    return group_fish, ratio

def heatmap(data, barcode):
    """
    This function inputs a DataFrame ('data') of reads and a barcode. It plots a heatmap showing
    the reads of that barcode across all fish, samples, and replicates.
    """
    data['sample'] = pd.to_numeric(data['sample'])
    group_fish, ratio = summarize_barcode(data, barcode)

    # sort fish in reverse order so that fish with large reads are on top.
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
