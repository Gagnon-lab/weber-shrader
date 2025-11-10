import numpy as np
import pandas as pd
import os
import ast
import json

import testis_model as tm

def check_scenario(scenario, overwrite = False):
    """
    This file checks whether the directory and data file for a given scenario already exists.

    If the directory does not exist, then one is created. If the file already exists,
    then an error is raised so that it is not overwritten. To suppress this error
    (meaning you are ok with overwriting data), set the optional parameter overwrite = True.
    """
    directory = f'./scenarios/{scenario}'
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not overwrite and os.path.exists(f'{directory}/{scenario}.csv'):
        raise Exception(f'A CSV file for scenario {scenario} already exists.')

def combine_and_save(data, scenario):
    """
    Takes a DataFrame `data` of results in tuple format and saves it to a file.
    

    """
    df = pd.concat(data, ignore_index = True, sort = False)

    for col in ['obs_result', 'hid_result', 'full_result']:
        if col in df.columns:
            df[col] = df[col].apply(lambda tup: tuple(int(x) for x in tup))
            df[col] = df[col].apply(str)
    
    df.to_csv(f'./scenarios/{scenario}/{scenario}.csv', index = False)

    return df

def find_scenarios():
    """
    This function searches the './scenarios' subdirectory and returns a
    list containing the names of all files that are found.
    """
    scenarios = []

    for scenario in os.listdir('scenarios'):
        path = os.path.join('scenarios', scenario, f'{scenario}.csv')
        if os.path.isfile(path):
            scenarios.append(scenario)

    return scenarios

def load_all_scenarios(scenarios):
    """
    Load data about fish and clones from a saved file.

    This function inputs the variable `group`, which represents the set of data to be loaded.
    The data loaded will come from './simulated-data/{group}.pkl').

    This function outputs two dictionaries: `clone_data` and `fish_data`. [TODO finish this]
    """
    clone_data = dict()
    fish_data = dict()

    if isinstance(scenarios, str):
        scenarios = [scenarios]
        
    for scenario in scenarios:
        df = pd.read_csv(f'./scenarios/{scenario}/{scenario}.csv')

        # Convert 'obs_result' from string to tuple of ints
        for col in ['obs_result', 'hid_result', 'full_result']:
            if col in df.columns:
                df[col] = df[col].apply(ast.literal_eval)
            
        fish_list = df['fish'].unique()
    
        for i, fish in enumerate(fish_list):
            fish_df = df.query('fish == @fish').copy().reset_index(drop = True)
            n_clones = len(fish_df['obs_result'].iloc[0])
    
            fish_data[scenario, fish] = fish_df
        
            for clone in range(n_clones):
                collapsed = tm.data_processing.collapse(fish_df, clone)
                collapsed.insert(2, 'clone', clone)
                clone_data[scenario, int(fish), clone] = collapsed

    return fish_data, clone_data

def load_clone_names():
    """
    Creates a dictionary `barcodes` such that `barcodes[fish, clone_name]`
    returns the index of that clone in the fish, while
    `barcodes[fish, index]` returns the name of the corresponding clone.
    """
    barcodes = dict()
    with open("./scenarios/real/barcodes.json", "r") as f:
        bc_file = json.load(f)
    
        for fish, bc_list in bc_file.items():
            for i, clone in enumerate(bc_list):
                barcodes[int(fish), i] = clone
                barcodes[int(fish), clone] = i