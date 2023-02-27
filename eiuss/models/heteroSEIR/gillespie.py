import numpy as np
from numba import jit
##################################################################
# model specific values
compartments = ['S', 'E', 'I_l', 'I_h', 'R']
infection_event = [['I_l', 'E', lambda x: x['beta'] * x['f_low']],
                   ['I_h', 'E', lambda x: x['beta'] * x['f_high']]]
transition_event = [['E', 'I_l', lambda x: x['gammaE']],  ## E1 -> I_l, I_h
                    ['E', 'I_h', lambda x: 0],
                    ['I_l', 'R', lambda x: x['gammaI']],
                    ['I_h', 'R', lambda x: x['gammaI']]]
##################################################################
compartments = {compartments[i]:i-1 for i in range(len(compartments))}      ## convert compartment to int index; idx of S is -1 as it is not included in genotype_count table

infection_src_idx  = np.array([compartments[x[0]] for x in infection_event])
infection_dst_idx  = np.array([compartments[x[1]] for x in infection_event])
transition_src_idx = np.array([compartments[x[0]] for x in transition_event])
transition_dst_idx = np.array([compartments[x[1]] for x in transition_event])



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
    statevar_next[1] -= np.sum(count_infection)
    _sum_by_index(statevar_next, count_infection, i_dst)

    # transition events
    _sum_by_index(statevar_next, (-1) * count_noninfection, n_src)
    _sum_by_index(statevar_next, count_noninfection[:-2],   n_dst)


    return




def transition_split(event_counts, idx_src, idx_dst, rate_dst, rng):
    event_counts[idx_dst] = rng.binomial(p=rate_dst, n=event_counts[idx_src])
    event_counts[idx_src] -= event_counts[idx_dst]

    return event_counts







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
        count = rng.poisson(rate_infection * dt * infectious_curr * S_curr)
        if np.sum(count) <= S_curr:
            break
        print(f'resimulating Gillespie step; may want to consider smaller dt: np.sum(count)({np.sum(count)}) > S_curr({S_curr})')
        print (f'rate = {rate_infection} * {dt} * {infectious_curr} * {S_curr}')

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