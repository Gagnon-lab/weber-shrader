import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.join(current_dir, '..')
sys.path.append(parent_dir)

import pipeline_utility


def get_anomaly_metrics(experiment, data_folder, root):
    """
    This function gets metrics for anomaly detection from data directory
    Metrics are returned as a touple where the first entry corresponds to GESTALT used entries and the second corresponds to reads remaining after filtering
    """
    info_df = pd.DataFrame({'sample' : [], 'usedEntries' : [], 'totalEntries' : [], 'dir' : []})
    # Concatenate merged.info files
    for dirpath, dirnames, filenames in os.walk(f'{root}/{data_folder}/mergedinfos'):
        for filename in filenames:          
            if filename.endswith('merged.info'):
                filepath = os.path.join(dirpath, filename)
                new_df = pd.read_csv(filepath, sep='\t')
                new_df['dir'] = dirpath
                info_df = pd.concat([info_df, new_df])
    total_counts = {}

    # Get total reads from .allReadCounts files
    for dirpath, dirnames, filenames in os.walk(f'{root}/{data_folder}/full-barcoded-fish'):
        for filename in filenames:
            if filename.endswith('.allReadCounts'):
                filepath = os.path.join(dirpath, filename)
                sample = filename.replace('.allReadCounts', '')
                read_counts = pd.read_csv(filepath, sep='\t')
                total_counts[sample] = read_counts['count'].sum()
    info_df['trueUsedEntries'] = info_df['sample'].apply(lambda x: total_counts.get(x, 0))

    # Create dataframe with GESTALT used entry information
    info_df['file'] = info_df['sample']
    info_df['usedEntryProp'] = info_df['usedEntries'] / info_df['totalEntries']
    info_df['trueUsedEntryProp'] = info_df['trueUsedEntries'] / info_df['totalEntries']
    info_df['fish'] = info_df['file'].apply(lambda x: x.split('_')[0])
    info_df['sample'] = info_df['file'].apply(lambda x: x.split('_')[1][1:])
    info_df['replicate'] = info_df['file'].apply(lambda x: x.split('_')[2][3:])
    info_df['sample'] = info_df['sample'].astype(int)

    # Get kept barcode info
    kept_file = f'{root}/experiments/{experiment}/prop_reads_kept.csv'
    kept_df = pd.read_csv(kept_file)
    kept_df.rename(columns={'proportion': 'keptProp'}, inplace=True)
    kept_df['sample'] = kept_df['sample'].astype(int)
    if not 'pass' in kept_df.columns:
        kept_df['pass'] = kept_df['pass_prop'] & kept_df['pass_z']
    return (info_df, kept_df)

def plot_fish(fish, ax, all_data):
    """
    Plots muller for a given fish
    """
    fish_data = all_data.query('fish == @fish').copy()
    pipeline_utility.plot_muller(all_data, fish, ax, map_samples=True)
    ax.set_title(f'{fish}')
    plt.tight_layout()


def show_anomalies(fish, axis, all_data, metrics, show_used=True, show_filtered=True, map_samples=False, show_cutoff=True):
    """
    Shows visualization for anomaly metrics for a given fish
    show_used: Shows GESTALT used entry proportion with red squares
    show_filered: Shows proportion of reads remaining after filtering pipeline with purple triangles
    map_samples: Maps the sample # to the month of the fish's life the sample was taken
    show_cutoff: Uses a dotted line to visualize the cutoff for a sample to be removed for each metric
    """
    fish_data = all_data.query('fish == @fish').copy()
    
    #Set colors
    used_entries_color = 'red'
    kept_values_color = 'purple'
    
    # Plot proportion of used and filtered entries
    if show_used:
        show_used_entries(metrics[0], fish, axis, used_entries_color, map_samples, show_cutoff)
    if show_filtered:
        show_filtered_proportion(metrics[1], fish, axis, kept_values_color, map_samples, show_cutoff)

    axis.set_xlabel('Time')
    axis.set_ylabel(None)
    xticks = np.unique(np.append(metrics[0][metrics[0]['fish'] == fish]['sample'], metrics[1][metrics[1]['fish'] == fish]['sample']))
    if map_samples:
        axis.set_xticks(np.array([pipeline_utility.sample_map[x] for x in xticks]))
    else:
        axis.set_xticks(xticks)
    axis.set_ylim(0, 1)

def show_used_entries(info_df, fish, axis, col, map_samples, show_cutoff):
    """
    Plots the proportion of reads used by GESTALT
    """
    fish_data = info_df.query('fish == @fish').copy()
    if map_samples:
        fish_data['sample'] = fish_data['sample'].map(lambda x: pipeline_utility.sample_map[x])
            
    fish_data['outlier'] = fish_data['trueUsedEntryProp'] < 0.3

    axis.plot(fish_data['sample'], fish_data['trueUsedEntryProp'], color=col, marker='s', linestyle='', markersize=5, alpha=0.5)
    axis.plot(fish_data[fish_data['outlier']]['sample'], fish_data[fish_data['outlier']]['trueUsedEntryProp'], color=col, marker='x', linestyle='', markersize=10)
    if show_cutoff:
        outlier = 0.3
        domain = np.array(range(fish_data['sample'].min(), fish_data['sample'].max() + 1))
        axis.plot(domain, np.repeat(outlier, len(domain)), color=col, linestyle='--', alpha=0.5)

def show_filtered_proportion(kept_df, fish, axis, col, map_samples, show_cutoff):
    """
    Plots the proportion of reads that remained after the filtering pipeline
    """
    fish_data = kept_df.query('fish == @fish').copy()
    if map_samples:
        fish_data['sample'] = fish_data['sample'].map(lambda x: pipeline_utility.sample_map[x])
    
    axis.plot(fish_data['sample'], fish_data['keptProp'], color=col, marker='^', linestyle='', markersize=5, alpha=0.5)
    axis.plot(fish_data[~fish_data['pass']]['sample'], fish_data[~fish_data['pass']]['keptProp'], color=col, marker='x', linestyle='', markersize=10)
    if show_cutoff:
        mean = fish_data['mean'].head(1).item()
        std = fish_data['std'].head(1).item()
        outlier = max(mean - 3 * std, 0.3)
        domain = np.array(range(fish_data['sample'].min(), fish_data['sample'].max() + 1))
        axis.plot(domain, np.repeat(outlier, len(domain)), color=col, linestyle='--', alpha=0.5)