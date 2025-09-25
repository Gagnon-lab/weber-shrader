import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.patches import Patch

def where_unedited(barcode_list):
    """
    Given a list of barcodes (strings), returns the index where the unedited barcode occurs.
    If the unedited barcode does not appear, return `None` instead.
    """
    unedited_barcode = 'NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE'
    where = np.where(barcode_list == unedited_barcode)[0]
    return where[0] if len(where) == 1 else None

def plot(data, fish, ax = None, replicate = None, add_legend = False, map_samples = False):
    sample_map = {1: 5, 2: 6, 3: 7, 4: 8, 5: 9, 6: 10, 7: 11, 8:12, 9:13, 10:14, 11:16, 12:17, 13:18, 14: 20, 15: 22, 16: 24, 17: 25, 18:26}
    subset = data.query('fish == @fish')
    if replicate is not None:
        subset = subset.query('replicate == @replicate')

    sorted_data = subset.sort_values(by=['sample', 'reads'], ascending=[True, False])
    barcode_order = sorted_data['barcode'].unique()
    
    grouped = subset.groupby(['barcode', 'sample'], as_index = False).aggregate({'reads': 'sum'})
    matrix = grouped.pivot(index = 'barcode', columns = 'sample', values = 'reads')
    matrix.fillna(0, inplace = True)
    matrix = matrix.loc[barcode_order]
    freqs = matrix.div(matrix.sum(axis=0), axis=1)
    
    # Source: https://github.com/BlakeRMills/MoMAColors/blob/main/R/Colors.R
    colors = np.array(["#ff4d6f", "#579ea4", "#df7713", "#f9c000", "#86ad34",
          "#5d7298", "#81b28d", "#7e1a2f", "#2d2651", "#c8350d", "#bd777a"])
    
    n_clones = np.shape(freqs)[0]
    repeats = -(-n_clones // len(colors))
    colors_extended = np.tile(colors, repeats)[:n_clones]
    
    unedited_index = where_unedited(freqs.index)
    if unedited_index is not None:
        # put in grey in the element corresponding to the unedited array.
        colors_extended = np.insert(colors_extended, unedited_index, '#aaaaaa')
    
        # remove the last index so that the array maintains the correct length.
        colors_extended = colors_extended[:-1]
    
    ts = matrix.columns
    if map_samples:
        ts = [sample_map[t] for t in ts]

    if ax is None:
        fig, ax = plt.subplots()
        
    ax.stackplot(ts, freqs[::-1], colors = colors_extended[::-1]);
    ax.set_xticks(ts);
    ax.set_xlabel("Time")
    ax.set_ylabel("Clonal frequency")
    ax.set_xlim(ts[0], ts[-1])
    ax.set_ylim(0, 1);

    # Create legend handles
    legend_labels = barcode_order
    legend_colors = colors_extended
    handles = [Patch(facecolor=c, label=b) for c, b in zip(legend_colors, legend_labels)]
    
    # Add legend to plot
    if add_legend:
        ax.legend(handles=handles[:20], title="Barcode", bbox_to_anchor=(1.05, 1), loc='upper left')

    if ax is None:
        return fig, ax