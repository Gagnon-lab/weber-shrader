import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

# Source: https://github.com/BlakeRMills/MoMAColors/blob/main/R/Colors.R
default_colors = np.array(["#ff4d6f", "#579ea4", "#df7713", "#f9c000", "#86ad34",
              "#5d7298", "#81b28d", "#7e1a2f", "#2d2651", "#c8350d", "#bd777a",
              '#1a5f7a', '#8b4c8f', '#4a9eff', '#ffd93d', '#3dcc8e', '#ff8c42',
              '#e63946', '#a0c334', '#6b4e71', '#00b4d8'])

def plot(results, ax = None, replicate = None, unedited_index = None, colors = None, offset_time = 0):
    """
    Input a `results` dataframe and `ax` axis and
    plot the corresponding Muller plot.

    If `replicate` is None, then all replicates are combined.
    Otherwise, `replicate` specifies which replicate will
    be plotted.
    """
    if ax is None:
        fig, ax = plt.subplots()

    if replicate is None:
        # For each sample, combine the results of all
        # replicates into a single row, where the
        # tuples representing reads are added
        # together elementwise.
        grouped = (
            results.groupby("time")["obs_result"]
            .agg(lambda tuples: tuple(map(sum, zip(*tuples))))
            .reset_index()
        )
    else:
        grouped = results[results.replicate == replicate]

    # obs_data has shape (n_clones, n_times)
    # Element i, j is the number of reads of clone i at time j.
    obs_data = np.array([list(result) for result in grouped.obs_result]).T

    # Convert reads to frequencies by dividing by row sums.
    freq_data = obs_data / np.sum(obs_data, axis = 0)

    # Extract all the times at which reads occur.
    times = np.array(grouped.time)

    if colors is None:
        colors = default_colors
        
    n_clones = np.shape(freq_data)[0]
    repeats = -(-n_clones // len(colors)) # same as n_clones/len(colors) rounded up
    colors_extended = np.tile(colors, repeats)[:n_clones]
    
    if unedited_index is not None:
        # put in grey in the element corresponding to the unedited array.
        colors_extended = np.insert(colors_extended, unedited_index, '#aaaaaa')

        # remove the last index so that the array maintains the correct length.
        colors_extended = colors_extended[:-1]

    times = times + offset_time
    times = np.round(times)
    ax.stackplot(times, *freq_data[::-1,:], colors=colors_extended[::-1])

    ax.set_xticks(times)
    
    ax.set_xlabel("Time")
    ax.set_ylabel("Barcode frequency")
    ax.set_xlim(times[0], times[-1])
    ax.set_ylim(0, 1);

def plot_tubules(spatial_result, time, ax=None, colors=None):
    """Plot the clonal composition of every tubule at a sample or record time."""
    if colors is None:
        colors = default_colors

    tubules = None
    for prefix in ("sample", "record"):
        times = np.asarray(spatial_result[f"{prefix}_times"])
        matches = np.flatnonzero(np.isclose(times, time))
        if len(matches):
            tubules = spatial_result[f"{prefix}_states"][matches[0]]
            break
    if tubules is None:
        available = np.unique(np.concatenate((spatial_result["sample_times"], spatial_result["record_times"])))
        raise ValueError(f"time {time} not found; available times are {available.tolist()}")
    if ax is None:
        _, ax = plt.subplots(figsize=(12, max(2.5, 0.6 * len(tubules))))

    max_clone = max(int(np.max(tubule)) for tubule in tubules if len(tubule))
    clone_colors = np.tile(colors, -(-(max_clone + 1) // len(colors)))
    row_spacing = 1.25
    for i, tubule in enumerate(tubules):
        y = i * row_spacing
        for x, clone in enumerate(tubule):
            ax.add_patch(Rectangle((x, y), 1, 1, facecolor=clone_colors[clone],
                                  ))#edgecolor="white", linewidth=0.25))

    max_length = max(map(len, tubules))
    ax.set(xlim=(0, max_length), ylim=(len(tubules) * row_spacing, 0),
           xlabel="Position", ylabel="Tubule", title=f"Time {time:g}")
    ax.set_yticks(np.arange(len(tubules)) * row_spacing + 0.5,
                  labels=np.arange(1, len(tubules) + 1))
    ax.set_aspect("auto")
    return ax
