import numpy as np
from numba import jit

##################################################################
# model specific values
compartments = ['S', 'E', 'I', 'R']
compartments = {compartments[i]:i-1 for i in range(len(compartments))}      ## convert compartment to int index; idx of S is -1 as it is not included in genotype_count table

infection_src_idx = np.array([compartments[x] for x in ['I']])
infection_dst_idx = np.array([compartments[x] for x in ['E']])

transition_src_idx = np.array([compartments[x] for x in ['E', 'I']])
transition_dst_idx = np.array([compartments[x] for x in ['I', 'R']])
##################################################################



# ------------------
# model dependent functions
@jit(nopython=True)
def update_statevar_next(statevar_curr, statevar_next, t_now_end, count_infection, count_noninfection,):

    i_dst = infection_dst_idx  + 2
    n_src = transition_src_idx + 2
    n_dst = transition_dst_idx + 2

    statevar_next[0] = t_now_end
    statevar_next[1:] = statevar_curr[1:]

    # infection events
    statevar_next[1] -= np.sum(count_infection)     ## susceptible
    #statevar_next[5] += np.sum(count_infection)
    _sum_by_index(statevar_next, count_infection, i_dst)


    # non-infection events
    _sum_by_index(statevar_next, (-1) * count_noninfection, n_src)
    _sum_by_index(statevar_next, count_noninfection[:-2],   n_dst)


    return


def assign_high_group (event, param_frac, rng, n_event_per_clade = 4, low_idx_ = 0, high_idx_ = 1):

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