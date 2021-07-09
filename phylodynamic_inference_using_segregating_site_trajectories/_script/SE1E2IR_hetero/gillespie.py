# ------------------
## gillespie.py (Salje_hetero; SE1E2IR model)
# ------------------
import numpy as np
from numba import jit
import sys

sys.path.append('..')
from tools import misc_ as misc

# ------------------
# model dependent functions
@jit(nopython=True)
def update_statevar_next(statevar_curr, statevar_next, t_now_end, count_infection, count_noninfection,):

    i_dst = np.array([0, 0, 0, 0]) + 2         # len(i_dst) = len(count_infection)         # todo -> put these in the param
    n_src = np.array([0, 1, 2, 0, 3, 4]) + 2   # len(n_src) = len(count_noninfection)
    n_dst = np.array([1, 2, 5, 3, 4, 5]) + 2   # len(n_dst) = len(count_noninfection) 
    """ """

    statevar_next[0] = t_now_end
    statevar_next[1:] = statevar_curr[1:]

    # infection events
    statevar_next[1] -= np.sum(count_infection)
    statevar_next[7] += np.sum(count_infection)
    _sum_by_index(statevar_next, count_infection, i_dst)


    # non-infection events
    _sum_by_index(statevar_next, (-1) * count_noninfection, n_src)
    _sum_by_index(statevar_next, count_noninfection[:-2], n_dst)


    return


def assign_high_group (event, param_frac, rng, n_event_per_clade = 6, low_idx_ = 0, high_idx_ = 3):

    i_event = np.arange(event.size)
    low_idx = i_event[i_event % n_event_per_clade == low_idx_]
    high_idx = i_event[i_event % n_event_per_clade == high_idx_]

    for i in range(len(low_idx)):
        l_idx = low_idx[i]; h_idx = high_idx[i]

        n_high = rng.binomial(p=param_frac, n=event[l_idx])

        event[l_idx] -= n_high
        event[h_idx] = n_high

    return event


# ------------------
# model independent functions
def event_count_noninfection (rate_noninfection, dt, statevar, rng):
    return _rng_poisson_vectorized(rate_noninfection, dt, statevar, statevar, rng)

def event_count_infection (rate_infection, dt, statevar_S, statevar, rng):
    S_curr = statevar_S
    infectious_curr = statevar

    if (S_curr < 0):
        print("ValueError: statevar < 0", (S_curr))

    if (infectious_curr < 0).all():
        print("ValueError: statevar < 0", (S_curr))

    while True:
        # l h /
        # E E / ....
        count = rng.poisson(rate_infection * dt * infectious_curr * S_curr)

        if np.sum(count) <= S_curr:
            break
        print('resimulating Gillespie step; may want to consider smaller dt')

    return count


def track_transmission_heterogeneity (transhet, count_infection, count_noninfection):

    # proportion of highly transmissible group
    E1_to_E2l = count_noninfection[0]
    E1_to_E2h = count_noninfection[3]

    # proportion of secondary infection caused by highly transmissible group - by E2 and I
    S_to_E_Il = count_infection[0] + count_infection[1]
    S_to_E_Ih = count_infection[2] + count_infection[3]


    transhet += np.array([1, E1_to_E2l, E1_to_E2h, S_to_E_Il, S_to_E_Ih ])

    return transhet




@jit(nopython=True)
def _sum_by_index(a, b, index):
    for i, ind in enumerate(index):
        a[ind] += b[i]
    return


def _rng_poisson(event_rate, dt, state_var, upper_limit, rng):

    if state_var < 0:
        print("ValueError: statevar < 0", (event_rate, dt, state_var))

    while True:
        k = rng.poisson(event_rate * dt * state_var)
        if k <= upper_limit:
            break
        print('resimulating Gillespie step; may want to consider decreasing dt')

    return k
_rng_poisson_vectorized = np.vectorize(_rng_poisson)