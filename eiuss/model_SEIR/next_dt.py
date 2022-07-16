import numpy as np

import utils
from model_SEIR import gillespie
from model_SEIR import genotype_count as genotype_count_
# ------------------
def update (particle, param, this_win_start, this_win_end):

    # clean up "recovered_dt"; remove individuals recovered during last sim_window
    particle.genotype_count[:,  param['idx_sampling']] = 0       

    # setup for for-loop
    t_idx = particle.last_idx       
    t_range = utils.np_arange_prevent_imprecision(this_win_end, this_win_start, param['dt'] * (-1))[::-1]
    # increase n_row of statevar if needed
    if (particle.last_idx + len(t_range) + 1) >= particle.statevar.shape[0]:
        particle.statevar = utils.add_row_to_nparray(particle.statevar, particle.last_idx + len(t_range))


    for idx in range(len(t_range)):
        t_now_end = t_range[idx]  # the end of current generation and the start of the next generation
        statevar_curr = particle.statevar[t_idx, :]
        statevar_next = particle.statevar[t_idx + 1, :]
        ## ----------------------------
        ## 1. get the number of events
        statevar_idx_i = [3]          # todo - multiclade
        statevar_idx_n = [2, 3]    # todo - multiclade

        event_count_infection  = gillespie.event_count_infection   (param['rate_infection'], param['dt'], statevar_curr[1], statevar_curr[statevar_idx_i], param['rng'])
        event_count_transition = gillespie.event_count_noninfection(param['rate_noninfection'], param['dt'], statevar_curr[statevar_idx_n], param['rng'])


        ## ----------------------------
        ## 2. update genotype_count
        genotype_count = particle.genotype_count
        n_genotype     = particle.n_genotype[0]
        
        ## 2-1. choose infection parents and mutate
        genotype_count_from_infection, n_new_genotype, particle.genotype_info, particle.n_site[0] = genotype_count_.infection_event_mutate (event_count_infection, statevar_curr, genotype_count, 0, particle.genotype_info, n_genotype, particle.n_site[0], param['mu_multinom'], param['rng'])
        
        ## 2-3. update genotype_count
        particle.genotype_count = genotype_count_.noninfection_event_update(event_count_transition, statevar_curr, particle.genotype_count, np.arange(0, 3), n_genotype, param['rng']);
        particle.genotype_count, particle.n_genotype[0] = genotype_count_.infection_event_update(genotype_count_from_infection.astype('int'), n_new_genotype.astype('int'), particle.genotype_count, np.arange(0, 3),  n_genotype)


        ## ----------------------------
        ## 3. update state_var
        statevar_next = particle.statevar[t_idx + 1, :]
        gillespie.update_statevar_next(statevar_curr, statevar_next, t_now_end, event_count_infection, event_count_transition)
        

        ## ----------------------------
        ## 4. safety check
        check_statevar_ = statevar_next[2:-1].astype("int")           ## except t, S, R
        check_genotype_ = particle.genotype_count[:, :-1].sum(axis = 0)  ## excetp R
        if np.array(check_statevar_ != check_genotype_).any():

            check_statevar = statevar_next.astype("int")           ## except t, S, R
            check_genotype = particle.genotype_count.sum(axis = 0)  ## excetp R

            error_str = f'------------------------------------------------\n' + \
                        f'>> check_statevar: {statevar_curr.astype("int")}\n' + \
                        f'>> check_statevar: {statevar_next.astype("int")}\n' + \
                        f'>> check_genotype: {particle.genotype_count.sum(axis = 0)}\n' + \
                        f'>> transition_eve: {event_count_transition}\n' + \
                        f'>> infection_even: {event_count_infection}\n' + \
                        f'------------------------------------------------\n' 

            raise Exception  (error_str)
            

        t_idx += 1

    particle.last_idx = t_idx

    return
