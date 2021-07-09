# ------------------
## next_dt.py (Salje_hetero; SE1E2IR model)
# ------------------
import numpy as np
import sys
sys.path.append('../')

from _script.tools import misc_ as misc

from _script.SE1E2IR_hetero import gillespie
from _script.SE1E2IR_hetero import genotype_count as genotype_count_
# ------------------
def update (particle, param, this_win_start, this_win_end):
    # parameters
    if 'mu_multinom' in param.keys():
        param_mu = param['mu_multinom']
    else:
        param_mu = param['mu']
    rate_infection = param['rate_infection']
    rate_noninfection = param['rate_noninfection']


    # idx - model dependent
    statevar_idx_i = np.array([1, 2, 3, 4]) + 2  # todo - multiclade
    statevar_idx_n = np.array([0, 1, 2, 0, 3, 4]) + 2  # todo - multiclade
    every_clade_3 = np.arange(param['n_clades']) * param['n_class'] + 5

    t_range = misc.np_arange_prevent_imprecision(this_win_end, this_win_start, param['dt'] * (-1))[::-1]
    # increase n_row of statevar if needed
    if (particle.last_idx + len(t_range) + 1) >= particle.statevar.shape[0]:
        particle.statevar = misc.add_row_to_nparray(particle.statevar, particle.last_idx + len(t_range))


    # clean up "recovered_dt"
    particle.genotype_count[:,  every_clade_3] = 0        # todo - any nicer way?


    t_idx = particle.last_idx
    for idx in range(len(t_range)):
        #print(t_idx)
        t_now_end = t_range[idx]  # the end of current generation and the start of the next generation
        statevar_curr = particle.statevar[t_idx, :]
        statevar_next = particle.statevar[t_idx + 1, :]

        ## get the number of transition events
        event_count_infection = gillespie.event_count_infection(rate_infection, param['dt'], statevar_curr[1], statevar_curr[statevar_idx_i], param['rng'])
        event_count_noninfection = gillespie.event_count_noninfection(rate_noninfection, param['dt'], statevar_curr[statevar_idx_n], param['rng'])
        event_count_noninfection = gillespie.assign_high_group(event_count_noninfection, param['p_H'], param['rng'])

        ## update genotype_count
        for i in range(param['n_clades']):
            clade_range = np.arange(i * param['n_class'], (i+1) * param['n_class'])
            genotype_count = particle.genotype_count[:, clade_range]
            n_genotype = particle.n_genotype[i]
            n_site = particle.n_site[i]

            genotype_count_from_infection, n_new_genotype, particle.genotype_info, particle.n_site[i] = genotype_count_.infection_event_mutate (event_count_infection, statevar_curr, genotype_count, i, particle.genotype_info, n_genotype, particle.n_site[i], param_mu, param['rng'])
            particle.genotype_count = genotype_count_.noninfection_event_update(event_count_noninfection, statevar_curr, particle.genotype_count, clade_range, n_genotype, param['rng']);
            particle.genotype_count, particle.n_genotype[i] = genotype_count_.infection_event_update(genotype_count_from_infection.astype('int'), n_new_genotype.astype('int'), particle.genotype_count, clade_range,  n_genotype)

        ## update state_var
        gillespie.update_statevar_next(statevar_curr, statevar_next, t_now_end, event_count_infection, event_count_noninfection)

        t_idx += 1


    particle.last_idx = t_idx

    return

def update_track_transhet (particle, param, this_win_start, this_win_end):
    # parameters
    if 'mu_multinom' in param.keys():
        param_mu = param['mu_multinom']
    else:
        param_mu = param['mu']
    rate_infection = param['rate_infection']
    rate_noninfection = param['rate_noninfection']


    # idx - model dependent
    statevar_idx_i = np.array([1, 2, 3, 4]) + 2  # todo - multiclade
    statevar_idx_n = np.array([0, 1, 2, 0, 3, 4]) + 2  # todo - multiclade
    every_clade_3 = np.arange(param['n_clades']) * param['n_class'] + 5

    t_range = misc.np_arange_prevent_imprecision(this_win_end, this_win_start, param['dt'] * (-1))[::-1]
    # increase n_row of statevar if needed
    if (particle.last_idx + len(t_range) + 1) >= particle.statevar.shape[0]:
        particle.statevar = misc.add_row_to_nparray(particle.statevar, particle.last_idx + len(t_range))


    # clean up "recovered_dt"
    particle.genotype_count[:,  every_clade_3] = 0        # todo - any nicer way?


    t_idx = particle.last_idx
    for idx in range(len(t_range)):
        #print(t_idx)
        t_now_end = t_range[idx]  # the end of current generation and the start of the next generation
        statevar_curr = particle.statevar[t_idx, :]
        statevar_next = particle.statevar[t_idx + 1, :]

        ## get the number of transition events
        event_count_infection = gillespie.event_count_infection(rate_infection, param['dt'], statevar_curr[1], statevar_curr[statevar_idx_i], param['rng'])
        event_count_noninfection = gillespie.event_count_noninfection(rate_noninfection, param['dt'], statevar_curr[statevar_idx_n], param['rng'])
        event_count_noninfection = gillespie.assign_high_group(event_count_noninfection, param['p_H'], param['rng'])

        ## update genotype_count
        for i in range(param['n_clades']):
            clade_range = np.arange(i * param['n_class'], (i+1) * param['n_class'])
            genotype_count = particle.genotype_count[:, clade_range]
            n_genotype = particle.n_genotype[i]
            n_site = particle.n_site[i]

            genotype_count_from_infection, n_new_genotype, particle.genotype_info, particle.n_site[i] = genotype_count_.infection_event_mutate (event_count_infection, statevar_curr, genotype_count, i, particle.genotype_info, n_genotype, particle.n_site[i], param_mu, param['rng'])
            particle.genotype_count = genotype_count_.noninfection_event_update(event_count_noninfection, statevar_curr, particle.genotype_count, clade_range, n_genotype, param['rng']);
            particle.genotype_count, particle.n_genotype[i] = genotype_count_.infection_event_update(genotype_count_from_infection.astype('int'), n_new_genotype.astype('int'), particle.genotype_count, clade_range,  n_genotype)

            # track transhet
            particle.transhet = gillespie.track_transmission_heterogeneity(particle.transhet, event_count_infection, event_count_noninfection)

        ## update state_var
        gillespie.update_statevar_next(statevar_curr, statevar_next, t_now_end, event_count_infection, event_count_noninfection)

        t_idx += 1


    particle.last_idx = t_idx

    return

