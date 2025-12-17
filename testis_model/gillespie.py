import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

from numba import njit

def verify_param_array(var, length):
    """
    This function inputs a variable `var` as well as its expected length
    (`length`) and its name (`name`).
    """
    if isinstance(var, (int, float)):  # Accept int or float
        var = np.full(length, var)
    var = np.asarray(var)
    if var.shape[0] != length:
        raise ValueError(f"Expected variable to have length {length}, but got {var.shape[0]}")
    return var

def sample_ic(n_stem, n_clones, rng = np.random.default_rng()):
    """
    Sample an initial condition for the model. This function assumes that there are
    `n_stem` stem cells and `n_clones` clones. The initial state is sampled from a
    uniform distribution, which can be derived by considering a simple birth
    process for each clone conditioned on a fixed total population.
    """
    if rng is None:
        rng = np.random.default_rng()
        
    samples = []

    # Choose n - 1 cut points from {1, ..., C - 1}
    cuts = rng.choice(range(1, n_stem), size=n_clones - 1, replace=False)
    cuts.sort()

    # Compute differences between cuts to get clone populations
    parts = [cuts[0]]
    for i in range(1, n_clones - 1):
        parts.append(cuts[i] - cuts[i - 1])
    parts.append(n_stem - cuts[-1])

    return np.array(parts)

class AdvancedTestisModel:
    def __init__(self, n_stem, n_clones, r_event, p_diff, r_decay, tau, sperm_per_diff, init_freq):
        """
        compute_rate: whether to compute the full rate transition matrix. This is necessary to run
        simulate().

        `init_freq`: initial condition of stem cell population. This can be input
        in two different ways:
        - To set exact clone populations at the start of the simulation,
          input `init_freq` as an array of length `n_clones` such that sum(init_freq)
          equals `n_stem`.
        - Alternative, `init_freq` can be input as an array of length `n_clones`
          such that sum(init_freq) equals 1. In this case, the initial clone
          populations will be set to a multinomial sample with probabilities
          given by `init_freq`.
        """
        self.n_stem = n_stem
        self.n_clones = n_clones
        self.init_freq = init_freq
        self.r_event = verify_param_array(r_event, n_clones)
        self.p_diff = verify_param_array(p_diff, n_clones)
        self.r_decay = r_decay
        self.sperm_per_diff = sperm_per_diff
        self.tau = tau

        self.diff_rates = self.r_event * self.p_diff

    def simulate_gillespie(self, sample_times, record_times=np.array([]),
                           init_time = 0, n_sperm=None, rng=None):
        if rng is None:
            rng = np.random.default_rng()
    
        state = np.zeros((2, self.n_clones), dtype=int)
        if np.sum(self.init_freq) == 1:
            state[0, :] = multinomial_numba(self.n_stem, self.init_freq, rng)
        elif np.sum(self.init_freq) == self.n_stem:
            state[0, :] = self.init_freq
        else:
            raise ValueError(f"init_freq must sum to 1 or n_stem, got {np.sum(self.init_freq)}")
    
        times = np.sort(np.concatenate((record_times, sample_times)))
        sample_flags = np.isin(times, sample_times)
        record_flags = np.isin(times, record_times)
    
        start_core = time.perf_counter()
    
        (
            sample_times_out,
            sample_stem,
            sample_sperm,
            sample_obs,
            sample_n_obs,
            record_times_out,
            record_stem,
            record_sperm,
            avg_stem_state
        ) = gillespie_core(
            init_time, times, state.copy(), self.r_event, self.p_diff,
            self.r_decay, self.tau, self.sperm_per_diff, self.diff_rates,
            self.n_clones, sample_flags, record_flags, n_sperm, rng
        )
    
        # Assemble sample_df
        sample_df = pd.DataFrame({
            'time': sample_times_out,
            'obs_result': [obs.tolist() for obs in sample_obs],
            'hid_result': [stem.tolist() for stem in sample_stem],
            'full_result': [sperm.tolist() for sperm in sample_sperm],
            'n_sperm': sample_n_obs
        })

        # Compute replicate
        replicate = sample_df.groupby("time").cumcount() + 1
        sample_df.insert(1, "replicate", replicate)
    
        # Assemble record_df
        record_df = pd.DataFrame({
            'time': record_times_out,
            'hid_result': [stem.tolist() for stem in record_stem],
            'full_result': [sperm.tolist() for sperm in record_sperm]
        })
    
        return sample_df, record_df, avg_stem_state

@njit
def weighted_choice(probs, rng):
    cdf = np.cumsum(probs)
    r = rng.random()
    return np.searchsorted(cdf, r)

@njit
def poisson_array(lam_array, rng):
    n = lam_array.shape[0]
    out = np.empty(n, dtype=np.int32)
    for i in range(n):
        out[i] = rng.poisson(lam_array[i])
    return out

@njit
def multinomial_numba(n, p, rng):
    out = np.zeros(len(p), dtype=np.int32)
    remaining = n
    total_p = np.sum(p)

    if total_p <= 0.0 or n <= 0:
        return out  # nothing to draw

    for i in range(len(p) - 1):
        if remaining <= 0:
            break
            
        pi = p[i] / np.sum(p[i:])
        pi = min(max(pi, 0.0), 1.0)    # clamp to [0, 1]
        draw = rng.binomial(remaining, pi)
        out[i] = draw
        remaining -= draw

    out[-1] = remaining
    return out

@njit
def dirichlet_numba(alpha_vec, rng):
    """
    Sample from Dirichlet(alpha_vec) using Numba.

    Parameters
    ----------
    alpha_vec : 1D array
        Positive concentration parameters (alpha * p).
    rng : np.random.Generator
        Numpy random number generator (must support .gamma).

    Returns
    -------
    sample : 1D array
        Probability vector sampled from Dirichlet(alpha_vec).
    """
    k = len(alpha_vec)
    gammas = np.empty(k)
    for i in range(k):
        gammas[i] = rng.gamma(shape=alpha_vec[i], scale=1.0)
    return gammas / gammas.sum()

@njit
def gillespie_core(init_time, times, state, r_event, p_diff, r_decay, tau, sperm_per_diff,
                   diff_rates, n_clones, sample_flags, record_flags, n_sperm, rng):

    t = init_time
    n_times = len(times)
    n_samples = np.sum(sample_flags)
    n_records = np.sum(record_flags)

    sample_times = np.empty(n_samples)
    sample_stem = np.empty((n_samples, n_clones), dtype = np.int32)
    sample_sperm = np.empty((n_samples, n_clones), dtype = np.int32)
    sample_obs = np.zeros((n_samples, n_clones), dtype = np.int32)
    sample_n_obs = np.zeros(n_samples, dtype = np.int32)
    # sample_obs = np.zeros((n_samples, n_clones), dtype = np.float64)

    record_times = np.empty(n_records)
    record_stem = np.empty((n_records, n_clones), dtype = np.int32)
    record_sperm = np.empty((n_records, n_clones), dtype = np.int32)

    sample_idx = 0
    record_idx = 0

    stem_accumulator = np.zeros(n_clones)
    total_time_accumulated = 0.0

    for i in range(n_times):
        target_time = times[i] - tau

        while t < target_time:
            stem_state = state[0, :]
            event_rate_per_clone = r_event * stem_state
            total_rate = np.sum(event_rate_per_clone)
            if total_rate == 0:
                break

            t_step = rng.exponential(1 / total_rate)
            if t_step >= target_time - t:
                break

            sperm_before = state[1, :].copy()

            j = weighted_choice(event_rate_per_clone / total_rate, rng)
            p_diff_j = p_diff[j]
            p_asymmetric = max(0.0, 1.0 - p_diff_j)
            outcome = weighted_choice(np.array([p_diff_j, p_asymmetric]), rng)

            if outcome == 0:
                div_probs = diff_rates * stem_state
                k = weighted_choice(div_probs / np.sum(div_probs), rng)
                state[0, j] -= 1
                state[1, j] += 2 * sperm_per_diff
                state[0, k] += 1
            # elif outcome == 1:
            #     diff_probs = diff_rates * stem_state
            #     k = weighted_choice(diff_probs / np.sum(diff_probs), rng)
            #     state[0, k] -= 1
            #     state[1, k] += 2 * sperm_per_diff
            #     state[0, j] += 1
            else:
                state[1, j] += sperm_per_diff

            lam_array = t_step * r_decay * sperm_before
            decayed_sperm = poisson_array(lam_array, rng)
            state[1, :] -= decayed_sperm
            t += t_step
            stem_accumulator += stem_state * t_step
            total_time_accumulated += t_step

        lam_array = (target_time - t) * r_decay * state[1, :]
        decayed_sperm = poisson_array(lam_array, rng)
        state[1, :] -= decayed_sperm
        t = target_time

        if record_flags[i]:
            record_times[record_idx] = target_time + tau
            record_stem[record_idx, :] = state[0, :]
            record_sperm[record_idx, :] = state[1, :]
            record_idx += 1

        if sample_flags[i]:
            sperm_total = np.sum(state[1, :])
            # if alpha is not None:
            if n_sperm is not None:
                # obs = state[1, :] / sperm_total
                if sperm_total == 0:
                    obs = np.zeros(n_clones, dtype = np.int32)
                else:
                    obs = multinomial_numba(n_sperm, state[1, :] / sperm_total, rng)
                    # param_vec = (alpha * state[1, :] / sperm_total).astype(np.float64)
                    # obs = dirichlet_numba(param_vec, rng)
                sample_obs[sample_idx, :] = obs
                sample_n_obs[sample_idx] = n_sperm
            else:
                sample_obs[sample_idx, :] = state[1, :]
                sample_n_obs[sample_idx] = sperm_total

            sample_times[sample_idx] = target_time + tau
            sample_stem[sample_idx, :] = state[0, :]
            sample_sperm[sample_idx, :] = state[1, :]
            sample_idx += 1

            # If next sample/record time is greater than current time,
            # remove all sperm.
            if i + 1 < n_times and times[i + 1] > times[i]:
                state[1, :] = 0

    avg_stem_state = stem_accumulator / total_time_accumulated

    return sample_times, sample_stem, sample_sperm, sample_obs, sample_n_obs, record_times, record_stem, record_sperm, avg_stem_state