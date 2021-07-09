# ------------------
## smc.py (Salje_hetero; SE1E2IR model)
# ------------------
import numpy as np
import sys


import statistics_S
import weight

def update_S(particle, n_samples, this_win_idx, R_idx, rng, count_for_mean = 50):

    for clade_idx in range(particle.n_clade):
        n_sample = n_samples[clade_idx]
        clade_R_idx = R_idx[clade_idx]

        segregating_site = []

        if n_sample < 1:    # if there is no sample, we can compute n_segregating
            particle.n_segregating[this_win_idx, clade_idx] = np.nan               # diff. from kk version

        else:
            n_genotype = particle.n_genotype[clade_idx]
            recovered_group = particle.genotype_count[:n_genotype, clade_R_idx]
            genotype_info = particle.genotype_info[:, 2 * clade_idx:2 * (clade_idx+1)]

            if np.sum(recovered_group) < n_sample:
                particle.n_segregating[this_win_idx, clade_idx] = np.nan

            else:
                for i in range(count_for_mean):
                    chosen_sample_idx = rng.choice(recovered_group.sum(), size=n_sample, replace=False)
                    chosen_sample_genotype = np.digitize(chosen_sample_idx, bins=recovered_group.cumsum(), right=False)

                    segregating_site.append(statistics_S.n_segregating_from_samples(chosen_sample_genotype, genotype_info))


                mean_segregating_site = sum(segregating_site) / count_for_mean
                particle.n_segregating[this_win_idx, clade_idx] = mean_segregating_site  # todo


        #return segregating_site
        return






def run_smc (particles, params, data, this_win_idx, this_win_start, this_win_end):

    this_w_vector = np.full([params['n_SMC_particles'], data['n_clades']], np.nan)
    no_introduction = False
    all_extinction = True

    data_flag = data['data_flag'][this_win_idx].astype('int')
    data_S = data['S_in_win'][this_win_idx]
    n_samples = (data['n_samples_in_win'][this_win_idx]).astype('int')

    # model dependent part
    cumI_idx = [x+2 for x in [5]]
    E1E2I_idx = [x+2 for x in [0, 1, 2, 3, 4]]
    R_idx = params['idx_g_R']


    for particle_idx in range(params['n_SMC_particles']):
        particle = particles[particle_idx]

        params['func_update_next_dt'](particle, params, this_win_start, this_win_end)       # update_next_dt
        update_S(particle, n_samples, this_win_idx, R_idx, params['rng'])      # todo

        t_idx = particle.last_idx
        # todo
        bool_cumI = particle.statevar[t_idx, cumI_idx] > 0
        bool_E_sum_I = [(particle.statevar[t_idx, 0+2] +\
                         particle.statevar[t_idx, 1+2] +\
                         particle.statevar[t_idx, 2+2] +\
                         particle.statevar[t_idx, 3+2] +\
                         particle.statevar[t_idx, 4+2]) > 0]


        for i in range(data['n_clades']):
            this_w_vector[particle_idx, i] = weight.get_weight_by_S(bool_cumI[i], bool_E_sum_I[i], data_flag[i],
                                                               particle.n_segregating[this_win_idx, i], data_S[i], np.nan)

        no_surviving_clade = (this_w_vector[particle_idx, :] == -1).all()
        all_extinction *= no_surviving_clade  # if any of particle has at least one surviving clade, all_extinction will be False

        if (this_w_vector[particle_idx,:] == -3).any():  # <time_intro inference>: abandon this time_intro as with this time_intro, the clade has not been introduced even when we have data
            no_introduction = True

    return this_w_vector, no_introduction, all_extinction





