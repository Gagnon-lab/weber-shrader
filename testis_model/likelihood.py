import numpy as np

from numba import njit
from numba.typed import List

from math import log, exp

from scipy.linalg import expm

@njit
def moran_rate_matrix_numba(n_stem, r):
    """
    This function computs the rate transition matrix for a Moran process
    with `n_stem` stem cells at rate r.
    """
    R = np.zeros((n_stem + 1, n_stem + 1))

    for i in range(n_stem + 1):
        rate_up = 0.0
        rate_down = 0.0

        if i > 0:
            rate_down = (n_stem - i) * i / n_stem
            R[i, i - 1] = rate_down
        if i < n_stem:
            rate_up = i * (n_stem - i) / n_stem
            R[i, i + 1] = rate_up

        R[i, i] = - (rate_up + rate_down)

    return r * R

@njit
def log_binom_pmf(k, n, p):
    if p == 0.0:
        return 0.0 if k == 0 else -np.inf
    if p == 1.0:
        return 0.0 if k == n else -np.inf

    # log binomial coefficient using log-factorials
    log_coeff = 0.0
    for i in range(1, k + 1):
        log_coeff += log(n - i + 1) - log(i)

    return log_coeff + k * log(p) + (n - k) * log(1 - p)

@njit
def compute_emission_matrix(n_stem, n_sperm):
    E = np.zeros((n_stem + 1, n_sperm + 1))
    for i in range(0, n_stem + 1):  # avoid p=0
        p = i / n_stem
        for j in range(n_sperm + 1):
            log_pmf = log_binom_pmf(j, n_sperm, p)
            E[i, j] = exp(log_pmf) if log_pmf > -np.inf else 0.0
    return E

#
# Add more variance when r is big
#

# from numba import njit
# import numpy as np
# from math import exp, log
# from math import lgamma

# @njit
# def log_beta_binom_pmf(k, n, a, b):
#     return (
#         lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
#         + lgamma(k + a) + lgamma(n - k + b)
#         - lgamma(n + a + b)
#         + lgamma(a + b) - lgamma(a) - lgamma(b)
#     )

# @njit
# def log_binom_pmf(k, n, p):
#     if p == 0.0:
#         return 0.0 if k == 0 else -np.inf
#     if p == 1.0:
#         return 0.0 if k == n else -np.inf
#     return (
#         lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
#         + k * log(p) + (n - k) * log(1 - p)
#     )

# @njit
# def compute_emission_matrix(n_stem, n_sperm, r):
#     E = np.zeros((n_stem + 1, n_sperm + 1))
#     for i in range(n_stem + 1):
#         p = i / n_stem
#         var = r * p * (1 - p) / n_stem
#         if var == 0.0:
#             # Use standard binomial
#             for j in range(n_sperm + 1):
#                 log_pmf = log_binom_pmf(j, n_sperm, p)
#                 E[i, j] = exp(log_pmf) if log_pmf > -np.inf else 0.0
#         else:
#             common = p * (1 - p) / var - 1
#             a = p * common
#             b = (1 - p) * common
#             for j in range(n_sperm + 1):
#                 log_pmf = log_beta_binom_pmf(j, n_sperm, a, b)
#                 E[i, j] = exp(log_pmf) if log_pmf > -np.inf else 0.0
#     return E

# @njit
# def forward_algorithm(init_state, T, E, obs):
#     n_states = T.shape[0]
#     n_obs = len(obs)

#     # Forward matrix: alpha[t, i] = P(obs_0:t, state_t = i)
#     alpha = np.zeros((n_obs, n_states))

#     # Initialization
#     for i in range(n_states):
#         alpha[0, i] = init_state[i] * E[i, obs[0]]

#     # Recursion
#     for t in range(1, n_obs):
#         for j in range(n_states):
#             total = 0.0
#             for i in range(n_states):
#                 total += alpha[t - 1, i] * T[i, j]
#             alpha[t, j] = total * E[j, obs[t]]

#     return alpha

@njit
def forward_algorithm(init_state, obs, Ts, Es, f_cutoff):
    n_states = Ts[0].shape[0]
    n_obs = len(obs)

    # Forward matrix: alpha[t, i] = P(obs_0:t, state_t = i)
    alpha = np.zeros((n_obs, n_states))
    log_scaling = np.zeros(n_obs)

    # Initialization
    for i in range(n_states):
        alpha[0, i] = init_state[i] * Es[0][i, obs[0]]

    # Scale alpha[0, :]
    max_val = np.max(alpha[0, :])
    if max_val > 0:
        alpha[0, :] /= max_val
        log_scaling[0] = np.log(max_val)
    else:
        log_scaling[0] = -np.inf  # or 0.0 if preferred
    
    # Recursion
    for t in range(1, n_obs):
        max_f = np.max(alpha[t - 1, :])
        threshold = max_f * f_cutoff
        active_rows = np.where(alpha[t - 1, :] >= threshold)[0]
        
        T = Ts[t]
        E = Es[t]
        for j in range(n_states):
            total = 0.0
            for i in active_rows:
                total += alpha[t - 1, i] * T[i, j]
            alpha[t, j] = total * E[j, obs[t]]

        # Scale alpha[t, :]
        max_val = np.max(alpha[t, :])
        if max_val > 0:
            alpha[t, :] /= max_val
            log_scaling[t] = np.log(max_val)
        else:
            log_scaling[t] = -np.inf

    return alpha, log_scaling

import numpy as np
from scipy.special import gammaln

def binomial_log_likelihood(p, samps, n_sperm):
    if p <= 0 or p >= 1:
        return -np.inf  # log-likelihood undefined at boundaries

    # Ensure inputs are arrays
    samps = np.asarray(samps)
    n_sperm = np.asarray(n_sperm)

    # Stable log binomial coefficient
    log_coeff = gammaln(n_sperm + 1) - gammaln(samps + 1) - gammaln(n_sperm - samps + 1)
    log_likelihood = log_coeff + samps * np.log(p) + (n_sperm - samps) * np.log(1 - p)

    return np.sum(log_likelihood)

from scipy.linalg import eig, inv
from scipy.stats import binom
from scipy.optimize import minimize, root_scalar

class LikelihoodModel:

    def __init__(self, n_stem):
        self.n_stem = n_stem
        self.T_cache = {}
        self.E_cache = {}
        self.R = moran_rate_matrix_numba(n_stem, 1)

        evals, evecs = eig(self.R)
        self.evals = np.real(evals)
        self.evecs = evecs
        self.evecs_inv = inv(evecs)

    def precompute_T(self, rs):
        for r in rs:
            self.get_T(r, save = True)

        self.rs = rs
   
    def get_T(self, r, save = False):
        # print(r)
        T = self.T_cache.get(r)
        if T is None:
            if r == 0:
                T = np.eye(self.n_stem+1)
            else:
                # T = expm(r * self.R)
                # Compute matrix exponentials using eigen decomposition
                scaled_vecs = self.evecs * np.exp(r * self.evals)
                T = scaled_vecs @ self.evecs_inv
            if save:
                self.T_cache[r] = T
        return T

    def get_E(self, n_sperm):
        E = self.E_cache.get(n_sperm)
        if E is None:
            E = compute_emission_matrix(self.n_stem, n_sperm)
            self.E_cache[n_sperm] = E
        return E

    def likelihood(self, samples, r, init_freq, f_cutoff = 1e-10, return_alpha = False, save = False,
                  init_time = 0):
        
        n_samples = len(samples)
        Ts = List()
        Es = List()

        if init_freq is None:
            init_freq = samples['obs_result'].iloc[0][0] / samples['n_sperm'].iloc[0]
            
        initial_state = binom.pmf(np.arange(self.n_stem + 1), n = self.n_stem, p = init_freq)

        times = samples['time']
        n_sperms = samples['n_sperm']

        t_gaps = [times[0] - init_time if i == 0 else times[i] - times[i-1] for i in range(len(times))]
        T_dict = {t_gap: self.get_T(t_gap * r, save = save) for t_gap in np.unique(t_gaps)}

        # print(np.unique(t_gaps))
        
        for i, t_gap in enumerate(t_gaps):
            Ts.append(T_dict[t_gap])
            Es.append(self.get_E(n_sperms[i]))

        obs = np.array(samples['obs_result'].apply(lambda x: x[0]))

        alpha, log_scaling = forward_algorithm(initial_state, obs, Ts, Es, f_cutoff)

        l = np.sum(alpha[-1,:])
        if not np.isfinite(l) or l <= 0:
            log_l = -np.inf
        else:
            log_l = np.sum(log_scaling) + np.log(l)
            
        if return_alpha:
            return log_l, alpha
        return log_l

    # def max_likelihood_null(self, samples):
    #     clone_samples = samples['obs_result'].apply(lambda x: x[0])
    #     n_sperms = samples['n_sperm']
    #     if_opt = clone_samples.sum() / n_sperms.sum()
        
    #     l_max = binomial_log_likelihood(if_opt, clone_samples, n_sperms)

    #     return if_opt, l_max

    def max_likelihood_null(self, samples, **kwargs):
        # this is numerically unstable
        clone_samples = samples['obs_result'].apply(lambda x: x[0])
        n_sperms = samples['n_sperm']
        ic = clone_samples.sum() / n_sperms.sum()

        def neg_likelihood(freq):
            ret_val = -self.likelihood(samples, 0, freq[0], **kwargs)
            # print(f'r = {r:.2f}, if = {init_freq:.2f}, like = {ret_val:.2f}.')
            
            return ret_val

        bounds = [(1e-2, 1 - 1e-2)]
        result = minimize(neg_likelihood, x0=ic, method='L-BFGS-B', bounds=bounds)
        
        if_opt  = result.x[0]
        l_max  = -result.fun

        return if_opt, l_max

    def max_likelihood(self, samples, statistic = 2, ic = None, **kwargs):
        # init_freq = samples['obs_result'].iloc[0][0] / samples['n_sperm'].iloc[0]
        # initial_state = binom.pmf(np.arange(n_stem + 1), n = lm.n_stem, p = init_freq)
        
        def neg_likelihood(params):
            r, init_freq = params
            ret_val = -self.likelihood(samples, r, init_freq, **kwargs)
            # print(f'r = {r:.2f}, if = {init_freq:.2f}, like = {ret_val:.2f}.')
            
            return ret_val

        if ic is None:
            ic = [36.00, 0.5]
            
        bounds = [(0.0, 1000.0), (1e-2, 1 - 1e-2)]
        
        result = minimize(neg_likelihood, x0=ic, method='L-BFGS-B', bounds=bounds)
        
        r_opt  = result.x[0]
        if_opt = result.x[1]
        l_max  = -result.fun
    
        if statistic is None:
            return r_opt, if_opt, l_max
            
        # Define target function for root finding
        def target(r):
            return self.likelihood(samples, r, if_opt, **kwargs) - (l_max - 2 * statistic)
    
        # Find lower bound root
        lower_root = None
        if r_opt > bounds[0][0] and target(0) < 0:
            try:
                lower_result = root_scalar(target, bracket=[bounds[0][0], r_opt], method='brentq')
                if lower_result.converged:
                    lower_root = lower_result.root
                else:
                    lower_root = 0.0
            except ValueError:
                lower_root = np.nan
        else:
            lower_root = 0.0
    
        # Find upper bound root
        upper_root = None
        if r_opt < bounds[0][1]:
            try:
                upper_result = root_scalar(target, bracket=[r_opt, bounds[0][1]], method='brentq')
                if upper_result.converged:
                    upper_root = upper_result.root
            except ValueError:
                pass
    
        return r_opt, if_opt, l_max, lower_root, upper_root
    
    def likelihood_grid(self, samples, rs, init_freq):
        ls = np.zeros_like(rs)

        for i, r in enumerate(rs):
            ls[i] = self.likelihood(samples, r, init_freq)

        return ls
        

    # def likelihood(self, data, r, initial_state):
    #     T = self.get_T(r)
    #     E = self.get_E(371)
    
    #     alpha = forward_algorithm(initial_state, T, E, data)
    #     return np.log(np.sum(alpha[-1,:]))