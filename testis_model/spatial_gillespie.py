import warnings

import numpy as np
import pandas as pd
from numba import njit

from testis_model.gillespie import multinomial_numba, poisson_array, verify_param_array, weighted_choice


def _choose_clone(remaining, candidates, mode, rng):
    if mode == "uniform":
        return int(rng.choice(candidates))
    weights = remaining[candidates]
    return int(rng.choice(candidates, p=weights / weights.sum()))


def _arrange_counts(counts, continuity, rng, clone_choice="uniform", first_clone=None):
    remaining = np.asarray(counts, dtype=np.int64).copy()
    arrangement = np.empty(remaining.sum(), dtype=np.int32)
    available = np.flatnonzero(remaining)
    if not len(available):
        return arrangement
    if first_clone is None or remaining[first_clone] == 0:
        current = _choose_clone(remaining, available, "remaining", rng)
    else:
        current = int(first_clone)

    for i in range(len(arrangement)):
        arrangement[i] = current
        remaining[current] -= 1
        if i == len(arrangement) - 1:
            break
        others = np.flatnonzero((remaining > 0) & (np.arange(len(remaining)) != current))
        if remaining[current] > 0 and (not len(others) or rng.random() < continuity):
            continue
        current = _choose_clone(remaining, others, clone_choice, rng)
    return arrangement


def initialize_spatial_state(counts, tubule_lengths, c1, c2, rng, clone_choice="uniform"):
    """Arrange fixed clone counts first globally with c1, then within tubules with c2."""
    global_state = _arrange_counts(counts, c1, rng, clone_choice, first_clone=0)
    spatial_state = np.empty_like(global_state)
    start = 0
    for length in tubule_lengths:
        stop = start + length
        local_counts = np.bincount(global_state[start:stop], minlength=len(counts))
        spatial_state[start:stop] = _arrange_counts(local_counts, c2, rng, clone_choice)
        start = stop
    return spatial_state


class SpatialAdvancedTestisModel:
    def __init__(self, n_stem, n_clones, n_tubules, r_event, p_diff, r_decay, tau,
                 sperm_per_diff, init_freq, c1, c2, clone_choice="uniform"):
        """
        Spatial neutral-drift model with linear, non-interacting tubules.

        ``init_freq`` may contain clone counts summing to ``n_stem`` or
        probabilities summing to one. ``clone_choice`` controls how a new clone
        is selected during initialization: ``"uniform"`` selects uniformly among
        eligible clone identities, while ``"remaining"`` weights by cells left.
        """
        if not isinstance(n_stem, (int, np.integer)) or not isinstance(n_tubules, (int, np.integer)):
            raise TypeError("n_stem and n_tubules must be integers")
        if n_tubules < 1 or n_stem < 2 * n_tubules:
            raise ValueError("n_stem must be at least 2 * n_tubules")
        if not 0 <= c1 <= 1 or not 0 <= c2 <= 1:
            raise ValueError("c1 and c2 must lie in [0, 1]")
        if clone_choice not in {"uniform", "remaining"}:
            raise ValueError("clone_choice must be 'uniform' or 'remaining'")

        self.n_stem = n_stem
        self.n_clones = n_clones
        self.n_tubules = n_tubules
        self.init_freq = np.asarray(init_freq)
        if self.init_freq.ndim != 1 or len(self.init_freq) != n_clones or np.any(self.init_freq < 0):
            raise ValueError("init_freq must be a nonnegative array of length n_clones")
        self.r_event = verify_param_array(r_event, n_clones).astype(float)
        self.p_diff = verify_param_array(p_diff, n_clones).astype(float)
        self.r_decay = r_decay
        self.sperm_per_diff = sperm_per_diff
        self.tau = tau
        self.c1 = c1
        self.c2 = c2
        self.clone_choice = clone_choice

        base, remainder = divmod(n_stem, n_tubules)
        self.tubule_lengths = np.full(n_tubules, base, dtype=np.int32)
        self.tubule_lengths[:remainder] += 1
        if remainder:
            warnings.warn(
                f"n_stem is not divisible by n_tubules; the first {remainder} tubule(s) "
                "receive one extra cell.", UserWarning, stacklevel=2
            )
        self.left_neighbor, self.right_neighbor = _neighbor_arrays(self.tubule_lengths)

    def _initial_counts(self, rng):
        total = self.init_freq.sum()
        if np.isclose(total, 1):
            return multinomial_numba(self.n_stem, self.init_freq.astype(float), rng)
        if np.isclose(total, self.n_stem) and np.allclose(self.init_freq, np.round(self.init_freq)):
            return self.init_freq.astype(np.int32)
        raise ValueError(f"init_freq must sum to 1 or n_stem, got {total}")

    def simulate_gillespie(self, sample_times, record_times=np.array([]), init_time=0,
                           n_sperm=None, rng=None):
        if rng is None:
            rng = np.random.default_rng()
        counts = self._initial_counts(rng)
        positions = initialize_spatial_state(
            counts, self.tubule_lengths, self.c1, self.c2, rng, self.clone_choice
        )
        sperm = np.zeros(self.n_clones, dtype=np.int64)
        times = np.sort(np.concatenate((record_times, sample_times)))
        sample_flags = np.isin(times, sample_times)
        record_flags = np.isin(times, record_times)

        result = spatial_gillespie_core(
            init_time, times, positions, counts, sperm, self.r_event, self.p_diff,
            self.r_decay, self.tau, self.sperm_per_diff, sample_flags, record_flags,
            n_sperm, self.left_neighbor, self.right_neighbor, rng
        )
        (sample_times_out, sample_stem, sample_sperm, sample_obs, sample_n_obs,
         sample_positions, record_times_out, record_stem, record_sperm,
         record_positions, avg_stem_state) = result

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
        spatial_result = {
            "tubule_lengths": self.tubule_lengths.copy(),
            "sample_times": sample_times_out.copy(),
            "sample_states": self._split_snapshots(sample_positions),
            "record_times": record_times_out.copy(),
            "record_states": self._split_snapshots(record_positions)
        }
        return sample_df, record_df, avg_stem_state, spatial_result

    def _split_snapshots(self, snapshots):
        cuts = np.cumsum(self.tubule_lengths)[:-1]
        return [[part.copy() for part in np.split(row, cuts)] for row in snapshots]


def _neighbor_arrays(tubule_lengths):
    n_stem = int(np.sum(tubule_lengths))
    left = np.full(n_stem, -1, dtype=np.int32)
    right = np.full(n_stem, -1, dtype=np.int32)
    start = 0
    for length in tubule_lengths:
        stop = start + int(length)
        left[start + 1:stop] = np.arange(start, stop - 1)
        right[start:stop - 1] = np.arange(start + 1, stop)
        start = stop
    return left, right


@njit
def spatial_gillespie_core(init_time, times, positions, stem_counts, sperm, r_event,
                           p_diff, r_decay, tau, sperm_per_diff, sample_flags,
                           record_flags, n_sperm, left_neighbor, right_neighbor, rng):
    t = init_time
    n_clones = len(stem_counts)
    n_samples = np.sum(sample_flags)
    n_records = np.sum(record_flags)
    sample_times = np.empty(n_samples)
    sample_stem = np.empty((n_samples, n_clones), dtype=np.int32)
    sample_sperm = np.empty((n_samples, n_clones), dtype=np.int64)
    sample_obs = np.zeros((n_samples, n_clones), dtype=np.int64)
    sample_n_obs = np.zeros(n_samples, dtype=np.int64)
    sample_positions = np.empty((n_samples, len(positions)), dtype=np.int32)
    record_times = np.empty(n_records)
    record_stem = np.empty((n_records, n_clones), dtype=np.int32)
    record_sperm = np.empty((n_records, n_clones), dtype=np.int64)
    record_positions = np.empty((n_records, len(positions)), dtype=np.int32)
    sample_idx = 0
    record_idx = 0
    stem_accumulator = np.zeros(n_clones)
    total_time_accumulated = 0.0

    for i in range(len(times)):
        target_time = times[i] - tau
        while t < target_time:
            event_rates = r_event[positions]
            total_rate = np.sum(event_rates)
            if total_rate <= 0:
                break
            t_step = rng.exponential(1 / total_rate)
            if t_step >= target_time - t:
                break

            stem_accumulator += stem_counts * t_step
            total_time_accumulated += t_step
            sperm -= poisson_array(t_step * r_decay * sperm, rng)
            t += t_step

            site = weighted_choice(event_rates / total_rate, rng)
            clone = positions[site]
            if rng.random() < p_diff[clone]:
                left = left_neighbor[site]
                right = right_neighbor[site]
                donor_site = right if left < 0 else left if right < 0 else (left if rng.random() < 0.5 else right)
                donor_clone = positions[donor_site]
                sperm[clone] += 2 * sperm_per_diff
                if donor_clone != clone:
                    positions[site] = donor_clone
                    stem_counts[clone] -= 1
                    stem_counts[donor_clone] += 1
            else:
                sperm[clone] += sperm_per_diff

        remaining_time = target_time - t
        if remaining_time > 0:
            stem_accumulator += stem_counts * remaining_time
            total_time_accumulated += remaining_time
            sperm -= poisson_array(remaining_time * r_decay * sperm, rng)
            t = target_time

        if record_flags[i]:
            record_times[record_idx] = target_time + tau
            record_stem[record_idx] = stem_counts
            record_sperm[record_idx] = sperm
            record_positions[record_idx] = positions
            record_idx += 1

        if sample_flags[i]:
            sperm_total = np.sum(sperm)
            if n_sperm is None:
                sample_obs[sample_idx] = sperm
                sample_n_obs[sample_idx] = sperm_total
            elif sperm_total == 0:
                sample_n_obs[sample_idx] = n_sperm
            else:
                sample_obs[sample_idx] = multinomial_numba(n_sperm, sperm / sperm_total, rng)
                sample_n_obs[sample_idx] = n_sperm
            sample_times[sample_idx] = target_time + tau
            sample_stem[sample_idx] = stem_counts
            sample_sperm[sample_idx] = sperm
            sample_positions[sample_idx] = positions
            sample_idx += 1
            if i + 1 < len(times) and times[i + 1] > times[i]:
                sperm[:] = 0

    avg_stem_state = stem_counts.astype(np.float64)
    if total_time_accumulated > 0:
        avg_stem_state = stem_accumulator / total_time_accumulated
    return (sample_times, sample_stem, sample_sperm, sample_obs, sample_n_obs,
            sample_positions, record_times, record_stem, record_sperm,
            record_positions, avg_stem_state)
