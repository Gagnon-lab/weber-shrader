import numpy as np
import pandas as pd
from numba import njit


def verify_param_array(var, length):
    if np.isscalar(var):
        var = np.full(length, var)
    var = np.asarray(var)
    if var.shape[0] != length:
        raise ValueError(f"Expected variable to have length {length}, but got {var.shape[0]}")
    return var


def sample_ic(n_stem, n_clones, rng=None, ka=0, kq=0):
    """Uniformly sample positive clone sizes, then split them into active and quiescent cells."""
    if rng is None:
        rng = np.random.default_rng()
    if n_clones < 1 or n_clones > n_stem:
        raise ValueError("n_clones must lie between 1 and n_stem")
    if ka < 0 or kq < 0:
        raise ValueError("ka and kq must be nonnegative")

    if n_clones == 1:
        clone_sizes = np.array([n_stem], dtype=np.int32)
    else:
        cuts = np.sort(rng.choice(np.arange(1, n_stem), size=n_clones - 1, replace=False))
        clone_sizes = np.diff(np.concatenate(([0], cuts, [n_stem]))).astype(np.int32)

    p_quiescent = kq / (ka + kq) if ka + kq > 0 else 0.0
    quiescent = np.zeros(n_clones, dtype=np.int32)
    if p_quiescent > 0:
        quiescent = rng.binomial(clone_sizes, p_quiescent).astype(np.int32)
    return np.vstack((clone_sizes - quiescent, quiescent))


class AdvancedTestisModel:
    def __init__(self, n_stem, n_clones, r_event, p_diff, r_decay, tau,
                 sperm_per_diff, init_freq, ka=0, kq=0):
        """
        Non-spatial stem-cell model with active and quiescent states.

        ``init_freq`` may be a one-dimensional clone probability/count vector or
        an exact array of shape ``(2, n_clones)``. Rows of the exact array are
        active and quiescent cells, respectively.
        """
        if ka < 0 or kq < 0:
            raise ValueError("ka and kq must be nonnegative")
        self.n_stem = n_stem
        self.n_clones = n_clones
        self.init_freq = np.asarray(init_freq)
        self.r_event = verify_param_array(r_event, n_clones)
        self.p_diff = verify_param_array(p_diff, n_clones)
        self.r_decay = r_decay
        self.sperm_per_diff = sperm_per_diff
        self.tau = tau
        self.ka = ka
        self.kq = kq
        self.diff_rates = self.r_event * self.p_diff

    def _initial_state(self, rng):
        if self.init_freq.shape == (2, self.n_clones):
            if np.any(self.init_freq < 0) or not np.allclose(self.init_freq, np.round(self.init_freq)):
                raise ValueError("A two-row init_freq must contain nonnegative integer counts")
            stem_state = self.init_freq.astype(np.int32)
            if stem_state.sum() != self.n_stem:
                raise ValueError(f"init_freq must contain {self.n_stem} total cells")
            return stem_state
        if self.init_freq.shape != (self.n_clones,) or np.any(self.init_freq < 0):
            raise ValueError("init_freq must have shape (n_clones,) or (2, n_clones)")

        total = self.init_freq.sum()
        if np.isclose(total, 1):
            clone_sizes = multinomial_numba(self.n_stem, self.init_freq.astype(float), rng)
        elif np.isclose(total, self.n_stem) and np.allclose(self.init_freq, np.round(self.init_freq)):
            clone_sizes = self.init_freq.astype(np.int32)
        else:
            raise ValueError(f"init_freq must sum to 1 or n_stem, got {total}")

        p_quiescent = self.kq / (self.ka + self.kq) if self.ka + self.kq > 0 else 0.0
        quiescent = np.zeros(self.n_clones, dtype=np.int32)
        if p_quiescent > 0:
            quiescent = rng.binomial(clone_sizes, p_quiescent).astype(np.int32)
        return np.vstack((clone_sizes - quiescent, quiescent))

    def simulate_gillespie(self, sample_times, record_times=np.array([]),
                           init_time=0, n_sperm=None, rng=None):
        if rng is None:
            rng = np.random.default_rng()
        stem_state = self._initial_state(rng)
        sperm_state = np.zeros(self.n_clones, dtype=np.int64)
        times = np.sort(np.concatenate((record_times, sample_times)))
        sample_flags = np.isin(times, sample_times)
        record_flags = np.isin(times, record_times)

        result = gillespie_core(
            init_time, times, stem_state, sperm_state, self.r_event, self.p_diff,
            self.r_decay, self.tau, self.sperm_per_diff, self.diff_rates,
            self.ka, self.kq, sample_flags, record_flags, n_sperm, rng
        )
        (sample_times_out, sample_stem, sample_sperm, sample_obs, sample_n_obs,
         record_times_out, record_stem, record_sperm, avg_stem_state) = result

        sample_df = pd.DataFrame({
            "time": sample_times_out,
            "obs_result": [x.tolist() for x in sample_obs],
            "hid_result": [x.tolist() for x in sample_stem],
            "full_result": [x.tolist() for x in sample_sperm],
            "n_sperm": sample_n_obs
        })
        sample_df.insert(1, "replicate", sample_df.groupby("time").cumcount() + 1)
        record_df = pd.DataFrame({
            "time": record_times_out,
            "hid_result": [x.tolist() for x in record_stem],
            "full_result": [x.tolist() for x in record_sperm]
        })
        return sample_df, record_df, avg_stem_state


@njit
def weighted_choice(probs, rng):
    return np.searchsorted(np.cumsum(probs), rng.random())


@njit
def poisson_array(lam_array, rng):
    out = np.empty(len(lam_array), dtype=np.int32)
    for i in range(len(lam_array)):
        out[i] = rng.poisson(lam_array[i])
    return out


@njit
def multinomial_numba(n, p, rng):
    out = np.zeros(len(p), dtype=np.int32)
    remaining = n
    if np.sum(p) <= 0 or n <= 0:
        return out
    for i in range(len(p) - 1):
        if remaining <= 0:
            break
        pi = min(max(p[i] / np.sum(p[i:]), 0.0), 1.0)
        out[i] = rng.binomial(remaining, pi)
        remaining -= out[i]
    out[-1] = remaining
    return out


@njit
def gillespie_core(init_time, times, stem_state, sperm_state, r_event, p_diff,
                   r_decay, tau, sperm_per_diff, diff_rates, ka, kq,
                   sample_flags, record_flags, n_sperm, rng):
    t = init_time
    n_clones = stem_state.shape[1]
    n_samples = np.sum(sample_flags)
    n_records = np.sum(record_flags)
    sample_times = np.empty(n_samples)
    sample_stem = np.empty((n_samples, 2, n_clones), dtype=np.int32)
    sample_sperm = np.empty((n_samples, n_clones), dtype=np.int64)
    sample_obs = np.zeros((n_samples, n_clones), dtype=np.int64)
    sample_n_obs = np.zeros(n_samples, dtype=np.int64)
    record_times = np.empty(n_records)
    record_stem = np.empty((n_records, 2, n_clones), dtype=np.int32)
    record_sperm = np.empty((n_records, n_clones), dtype=np.int64)
    sample_idx = 0
    record_idx = 0
    stem_accumulator = np.zeros((2, n_clones))
    total_time_accumulated = 0.0

    for i in range(len(times)):
        target_time = times[i] - tau
        while t < target_time:
            active = stem_state[0]
            quiescent = stem_state[1]
            event_rates = np.concatenate((r_event * active, kq * active, ka * quiescent))
            total_rate = np.sum(event_rates)
            if total_rate == 0:
                break
            t_step = rng.exponential(1 / total_rate)
            if t_step >= target_time - t:
                break

            stem_accumulator += stem_state * t_step
            total_time_accumulated += t_step
            sperm_before = sperm_state.copy()
            event = weighted_choice(event_rates / total_rate, rng)

            if event < n_clones:
                clone = event
                if rng.random() < p_diff[clone]:
                    division_rates = diff_rates * active
                    replacement = weighted_choice(division_rates / np.sum(division_rates), rng)
                    active[clone] -= 1
                    sperm_state[clone] += 2 * sperm_per_diff
                    active[replacement] += 1
                else:
                    sperm_state[clone] += sperm_per_diff
            elif event < 2 * n_clones:
                clone = event - n_clones
                active[clone] -= 1
                quiescent[clone] += 1
            else:
                clone = event - 2 * n_clones
                quiescent[clone] -= 1
                active[clone] += 1

            sperm_state -= poisson_array(t_step * r_decay * sperm_before, rng)
            t += t_step

        remaining_time = target_time - t
        if remaining_time > 0:
            stem_accumulator += stem_state * remaining_time
            total_time_accumulated += remaining_time
            sperm_state -= poisson_array(remaining_time * r_decay * sperm_state, rng)
            t = target_time

        if record_flags[i]:
            record_times[record_idx] = target_time + tau
            record_stem[record_idx] = stem_state
            record_sperm[record_idx] = sperm_state
            record_idx += 1

        if sample_flags[i]:
            sperm_total = np.sum(sperm_state)
            if n_sperm is None:
                sample_obs[sample_idx] = sperm_state
                sample_n_obs[sample_idx] = sperm_total
            elif sperm_total == 0:
                sample_n_obs[sample_idx] = n_sperm
            else:
                sample_obs[sample_idx] = multinomial_numba(n_sperm, sperm_state / sperm_total, rng)
                sample_n_obs[sample_idx] = n_sperm
            sample_times[sample_idx] = target_time + tau
            sample_stem[sample_idx] = stem_state
            sample_sperm[sample_idx] = sperm_state
            sample_idx += 1
            if i + 1 < len(times) and times[i + 1] > times[i]:
                sperm_state[:] = 0

    avg_stem_state = stem_state.astype(np.float64)
    if total_time_accumulated > 0:
        avg_stem_state = stem_accumulator / total_time_accumulated
    return (sample_times, sample_stem, sample_sperm, sample_obs, sample_n_obs,
            record_times, record_stem, record_sperm, avg_stem_state)
