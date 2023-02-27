import numpy as np
from numba import jit
##################################################################
# model specific values
compartments     = ['S', 'E1', 'E2_l', 'I_l', 'E2_h', 'I_h', 'R', 'O_l', 'O_h']
infection_event  = [['E2_l', 'E1', lambda x: x['beta'] * x['f_low' ]],
                    [ 'I_l', 'E1', lambda x: x['beta'] * x['f_low' ]],
                    ['E2_h', 'E1', lambda x: x['beta'] * x['f_high']],
                    [ 'I_h', 'E1', lambda x: x['beta'] * x['f_high']],
                    [ 'O_l', 'E1', lambda x: x['beta'] * x['f_low' ] * x['eta']],
                    [ 'O_h', 'E1', lambda x: x['beta'] * x['f_high'] * x['eta']]]

transition_event = [['E1'   , 'E2_l',   lambda x: x['gammaE1']],
                    ['E2_l' ,  'I_l',   lambda x: x['gammaE2']],
                    [ 'I_l' ,  'R'  ,   lambda x: x['gammaI' ]],
                    ['E1'   , 'E2_h',   lambda x: 0           ],        ## E1 -> E2_l and E2_h; binomial distribution
                    ['E2_h' ,  'I_h',   lambda x: x['gammaE2']],
                    [ 'I_h' ,  'R'  ,   lambda x: x['gammaI' ]]]
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
        if ((np.sum(count) <= S_curr) * (count <= infectious_curr)).all():
        #if (np.sum(count) <= S_curr):
            break
        print(f'resimulating Gillespie step; may want to consider smaller dt: np.sum(count)({np.sum(count)}) > S_curr({S_curr})')
        print (f'rate = {rate_infection} * {dt} * {infectious_curr} * {S_curr}')

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