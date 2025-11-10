import pandas as pd
from saving_loading import find_scenarios, load_all_scenarios
from joblib import Parallel, delayed
from tqdm import tqdm
import argparse
import os
from itertools import islice
import traceback
import time
import warnings

from metrics import init_time, freq_at, run_likelihood
import testis_model as tm

lm = tm.likelihood.LikelihoodModel(1000)

def process_clone(data):
    series = dict()
    it = init_time(data)
    if it is None:
        # There is no observed data from this clone.
        return None
    series['init_time'] = it
    series['init_freq'] = freq_at(data, series['init_time'])

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=RuntimeWarning,
            message=".*invalid value encountered in subtract.*"
        )
        r_opt, if_opt, l_max, lower, upper, if_null, l_null = \
            run_likelihood(data, lm = lm)

    series['r_opt'] = r_opt
    series['if_opt'] = if_opt
    series['l_max'] = l_max
    series['lower'] = lower
    series['upper'] = upper
    series['if_null'] = if_null
    series['l_null'] = l_null

    return series

def process_scenario(scenario, n_sims, overwrite):
    start_time = time.time()
    directory = f'./scenarios/{scenario}'

    if n_sims is None:
        filename = f'{directory}/{scenario}_long.csv'
    else:
        filename = f'{directory}/{scenario}_long_test.csv'

    if not os.path.isdir(directory):
        return {'scenario': scenario, 'status': 'error - missing directory'}

    if os.path.exists(filename) and not overwrite:
        return {'scenario': scenario, 'status': 'skipped - file already exists'}

    fish_data, clone_data = load_all_scenarios(scenario)

    if n_sims is None:
        items = clone_data.items()
    else:
        items = islice(clone_data.items(), n_sims)
    
    all_series = []
    for (scenario, fish, clone), data in items:
        try:    
            series = process_clone(data)
            if series is None:
                continue
            series['scenario'] = scenario
            series['fish'] = fish
            series['clone'] = clone
            
            all_series.append(series)
        except Exception as e:
            tb = traceback.format_exc()

            duration = time.time() - start_time
            return {'scenario': scenario, 'status': f'crash in fish {fish} clone {clone}, see traceback. Total runtime: {duration:.2f}s.', 'traceback': tb}
    
    df = pd.DataFrame(all_series)
    df.to_csv(filename, index = False)
    
    duration = time.time() - start_time
    return {'scenario': scenario, 'status': f'Success! Results saved to {filename}. Total runtime: {duration:.2f}s.'}

def main(min_idx, max_idx, n_sims, overwrite):
    print('Loading available scenarios...')
    scenarios = find_scenarios()
    print(f'There are {len(scenarios)} found scenarios: {scenarios}')
    max_idx = min(max_idx, len(scenarios))
    print(f'I will run scenarios[{min_idx}:{max_idx}].')

    scenarios_to_run = scenarios[min_idx:max_idx]

    results = Parallel(n_jobs=-1)(
        delayed(process_scenario)(s, n_sims, overwrite) for s in tqdm(scenarios_to_run, desc="Processing scenarios")
    )

    print('Finished running scenarios. Results are: ')
    for res in results:
        scenario = res['scenario']
        status = res['status']
        print(f'{scenario:<20}{status}')

        if 'traceback' in res:
            print(res['traceback'])
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run long_metrics with optional scenarios.')
    
    parser.add_argument(
        '-s', '--scenarios',
        nargs = 2,
        type=int,
        metavar=('MIN_IDX', 'MAX_IDX'),
        help = 'Range of scenario names to run [min, max)'
    )
    
    parser.add_argument(
        '-o', '--overwrite',
        action = 'store_true',
        help = 'Overwrite existing outputs (default: False)'
    )

    parser.add_argument(
        '-n', '--nsims',
        type=int,
        default=None,
        help='Number of simulations to run (default: None)'
    )
    
    args = parser.parse_args()
    min_idx, max_idx = args.scenarios
    main(min_idx, max_idx, args.nsims, args.overwrite)