import os
import pandas as pd
import numpy as np
import colorsys

def hue(r):
    """
    This function inputs a drift rate `r` and outputs a shade of blue.
    Lower values of `r` correspond to darker shades.

    This function is used to make colors consistent across the figures.
    """
    lamb = (symlog10(r) - symlog10(0)) / (symlog10(36) - symlog10(0))
    return colorsys.hls_to_rgb(0.5, 0.30 + 0.45 * lamb, 0.5)

def symlog10(x, linthresh=1e-2):
    return np.sign(x) * np.log10(1 + np.abs(x) / linthresh)

def inv_symlog10(y, linthresh=1e-2):
    return np.sign(y) * ((10**np.abs(y) - 1) * linthresh)

exp_hue = '#ffc566'
mp_hue = '#f75973'
# trunc_hue = '#c56eeb'
tq_hue = '#c56eeb'
no_collapse_hue = '#3ddb6f'

small_label_settings = {
    'font.family': 'Arial',
    'font.size': 7,              # base font size
    'xtick.labelsize': 6,        # x-axis tick label size
    'ytick.labelsize': 6,        # y-axis tick label size
    'axes.labelsize': 7,         # axis label size
    'axes.titlesize': 7,         # title size
    'legend.fontsize': 7,        # legend font size
    'figure.titlesize': 7,        # figure title size
}

large_label_settings = {
    'font.family': 'Arial',
    'font.size': 14,              # base font size
    'xtick.labelsize': 9,        # x-axis tick label size
    'ytick.labelsize': 9,        # y-axis tick label size
    'axes.labelsize': 12,         # axis label size
    'axes.titlesize': 12,         # title size
    'legend.fontsize': 12,        # legend font size
    'figure.titlesize': 12,        # figure title size
}

palette = {
    'Model\n($rp=0$)': hue(0),
    'Model ($rp=0$)': hue(0),
    'Model': hue(2),
    'Model\n($rp=2$)': hue(2),
    'Model ($rp=2$)': hue(2),
    'Model\n($rp=0.5$)': hue(0.5),
    'Model\n($rp=8$)': hue(8),
    'Model\n($rp=25$)': hue(25),
    'Exp. data': exp_hue,
    'Experimental data': exp_hue,
    'Exp. data\n(MP 0.95)': mp_hue,
    'Exp. data (MP 0.95)': mp_hue,
    'Exp. data\n(TQ 14)': tq_hue,
    'Exp. data (TQ 14)': tq_hue,
    'Exp. data\n(no cluster)': no_collapse_hue,
    'Exp. data (no cluster)': no_collapse_hue,
    # 'Exp. data\n(trunc.)': trunc_hue,
    # 'Exp. data (truncated)': trunc_hue
}

# master_fish_order and master_colors define the order and colors for many
# main and supplemental figure panels.
master_fish_order = [
    "fish1", "fish2", "fish4", "fish5", "fish8",
    "fish10", "fish12", "fish13", "fish14", "fish15",
    "fish16", "fish17", "fish18", "fish22", "fish26",
    "fish35", "fish37", "fish38", "fish40", 
    "fish42", "fish44", "fish45", "fish50", 
    "fish51" 
]

master_colors = [
    "#14466C", "#2278B5", "#add8e6", #blues
    "#CD6C18", "#F57F20", "#FCBB78", #oranges
    "#1A702E", "#2FA148", "#9ED08A", #greens
    "#C4AA00", "#FCEA00", "#FFF468", #yellows
    "#940C07","#E62923", "#F69797", #reds
    "#673484", "#9268AD", "#C5B1D6", #purples
    "#CD2A91", "#E580BD", "#FCC7DC", #pinks
    "#008080", "#43b3ae", "#96ded1", #teals
    "#282B3D", "#7C83BC", "#C3CDE6" #periwinkles
]

color_map = dict(zip(master_fish_order, master_colors))

# Load in the clone results file.

file_path = os.path.join('..', 'model_simulations', 'all_clone_results.csv')
df = pd.read_csv(file_path)

exp_scenario = 'full_pipeline'

# Only include data that has at least 2 observations and are not trivial.
invalid_real_entries = df.query('scenario == @exp_scenario').query('n_transient < 2')
df = df.query('n_transient > 1').copy()

df['l_gap'] = df['l_max'] - df['l_null']

# Test statistic 3.84 corresponds to
# an alpha=0.05 signifiance level.
df['pass_null'] = 2 * df['l_gap'] > 3.84

with np.errstate(divide='ignore', invalid='ignore'):
    df['log_r'] = np.log(df['r_opt'])
    df['log10_r'] = np.log10(df['r_opt'])
    
df['symlog10_r'] = symlog10(df['r_opt'])