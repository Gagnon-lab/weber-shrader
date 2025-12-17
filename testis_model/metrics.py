import numpy as np
import pandas as pd
import testis_model as tm

def freq(clone_obs):
    return clone_obs[0] / (clone_obs[0] + clone_obs[1])

def init_time(data):
    time = data.loc[data['obs_result'].apply(lambda x: x[0] != 0), 'time'].min()
    if np.isnan(time):
        return None
    else:
        return data.loc[data['obs_result'].apply(lambda x: x[0] != 0), 'time'].min()

# def final_time(data):
#     return int(data.loc[data['obs_result'].apply(lambda x: x[0] != 0), 'time'].max())

def min_time(data):
    return data['time'].min()

def max_time(data):
    return data['time'].max()

def last_time(data):
    """
    Returns the last time before the clone is lost (OR is fixed).
    """
    max_t = max_time(data)
    final_alive_t = data.loc[data['obs_result'].apply(lambda x: x[0] > 0 and x[1] > 0), 'time'].max()

    return final_alive_t

def n_transient(data):
    """
    Returns the number of time points not at fixation or loss.
    """
    return data.loc[data['obs_result'].apply(lambda x: x[0] > 0 and x[1] > 0), 'time'].nunique()

def freq_at(data, t):
    state_at_t = tm.data_processing.add_tuples_elementwise(
        data.query('time == @t')['obs_result'])
    return freq(state_at_t)

def size(data):
    total_state = tm.data_processing.add_tuples_elementwise(data.obs_result)
    return freq(total_state)

# def sum_tuples(series):
#     return tm.data_processing.add_tuples_elementwise(*series)

def avg_freq_change(data):
    grouped = data.groupby('time', as_index=False).agg({'obs_result': tm.data_processing.add_tuples_elementwise})
    grouped['freq'] = grouped['obs_result'].apply(freq)
    grouped.set_index('time', inplace = True)
    # print(grouped)
    
    it = init_time(data)
    dt = last_time(data)
    if np.isnan(dt):
        return np.nan

    if dt == it:
        return np.nan
    else:
        return (grouped.loc[it]['freq'] - grouped.loc[dt]['freq']) / (dt - it)

def run_likelihood(data, lm = None, statistic = 4):
    if lm is None:
        lm = tm.likelihood.LikelihoodModel(1000)

    init_t = init_time(data) - lm.tau

    # Only use data after the clone appears
    used_data = data[data['time'] >= tm.metrics.init_time(data)].reset_index()
        
    # data = tm.data_processing.normalize_sperm_counts(data.copy(), 371)
    r_opt, l_max, lower, upper = lm.max_likelihood(used_data, statistic = statistic, init_time = init_t)
    l_null = lm.likelihood(used_data, 0, None, init_time = init_t, f_cutoff = 0)

    return r_opt, l_max, lower, upper, l_null
        
# def run_binom_likelihood(data, test_statistic = 4):
#     data = tm.data_processing.normalize_sperm_counts(data.copy(), 500)
#     it = init_time(data)
#     init_size = freq_at(data, it)
#     init_sperm_freq = np.array([init_size, 1 - init_size])

#     def likelihood_binom(r):
#         params = dict(
#             n_stem = 200, n_clones = 2, lamb = r,
#             init_freq = init_sperm_freq, compute_rate = True,
#             compute_emission = True, dense = True
#         )
#         with np.errstate(divide='ignore', invalid='ignore'):
#             hmm = tm.moran.Moran(**params)
#             return np.log(hmm.likelihood(data))

#     with np.errstate(divide='ignore', invalid='ignore'):
#         r_binom, lb_binom, ub_binom, ml_binom = \
#             tm.max_likelihood.max_likelihood_and_bounds(likelihood_binom, test_statistic)
    
#     likelihood_0_binom = likelihood_binom(0)

#     return r_binom, lb_binom, ub_binom, ml_binom, likelihood_0_binom

# def run_beta_likelihood(data, test_statistic = 4):
#     data = tm.data_processing.normalize_sperm_counts(data.copy(), 500)
#     it = init_time(data)
#     init_size = freq_at(data, it)
#     init_sperm_freq = np.array([init_size, 1 - init_size])

#     def likelihood_beta(r):
#         params = dict(
#             n_stem = 200, n_clones = 2, lamb = r,
#             init_freq = init_sperm_freq, compute_rate = True,
#             compute_emission = True, dense = True
#         )
#         with np.errstate(divide='ignore', invalid='ignore'):
#             hmm = tm.moran.MoranBeta(**params)
#             return np.log(hmm.likelihood(data))

#     with np.errstate(divide='ignore', invalid='ignore'):
#         r_beta, lb_beta, ub_beta, ml_beta = \
#             tm.max_likelihood.max_likelihood_and_bounds(likelihood_beta, test_statistic)
    
#     likelihood_0_beta = 0#likelihood_beta(0)

#     return r_beta, lb_beta, ub_beta, ml_beta, likelihood_0_beta