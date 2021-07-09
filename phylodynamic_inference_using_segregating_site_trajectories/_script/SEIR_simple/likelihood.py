# ------------------
## likelihood.py (SEIR without tranmission heterogeneity)
# ------------------
import copy
import numpy as np
import sys

from SEIR_simple import smc
# ------------------
def calculate (particles, params, imported_data, config):

    params = config.get_infection_rate(params)

    # get simulation window
    sim_windows = np.arange(max(imported_data['win_intervals']), params['timeintro'], imported_data['window_dt'] * (-1))
    sim_windows = sim_windows[::-1]
    data = _match_data_and_simwin(imported_data, sim_windows, imported_data['n_clades'])



    # initialize particles & add n_segregating row if needed
    for particle in particles:
        particle.initiate(params, len(sim_windows))
        particle.n_segregating[:len(data['win_intervals']), -1] = data['win_intervals']


    # matrix for weight, resampling
    w_matrix = np.full([len(sim_windows), params['n_SMC_particles']], np.nan)                 # shape = sim_window * n_SMC_particles
    m_matrix = np.full([len(sim_windows), params['n_SMC_particles']], np.nan)                 # weight of "meaningful" particles (weight >= 0)
    k_matrix = np.full([len(sim_windows), params['n_SMC_particles']], -9999, 'int')           # shape  = sim_window * n_SMC_particles


    # simulate and resample particles for each window
    this_win_start = params['timeintro']
    for this_win_idx in range(len(sim_windows)):
        this_win_end = sim_windows[this_win_idx]
        print("loglikelihood :", (this_win_idx, this_win_start, this_win_end))

        this_w_vector, no_introduction, all_extinction = smc.run_smc(particles, params, data, this_win_idx, this_win_start, this_win_end)


        if no_introduction :
            # when infering time_intro,
            # if a clade is not introduced in any of the particle but there is a data exist,
            # this time_intro should be abandoned and the log_likelihood = - inf
            log_likelihood = -np.inf
            print ("no introduction")
            return log_likelihood, particles, True


        if all_extinction :
            print ("all_extinction")
            # if all particles went extinction in this window interval,
            # nothing will be changed in the next window interval
            # -> we want at least one particle with at least one surviving clade
            # -> rerun for this window interval until we get "at least one particle with at least one surviving clade"



        this_w_vector_prod = np.prod(this_w_vector, axis = 1)
        w_matrix[this_win_idx, :] = np.prod(this_w_vector, axis=1)
        k_matrix[this_win_idx, :] = _resample_particle(np.where(this_w_vector_prod < 0, 0, this_w_vector_prod), params['n_SMC_particles'])


        resample_unique, resample_count = np.unique(k_matrix[this_win_idx, :], return_counts=True)
        resample_repeated_idx_in_unique = np.repeat(np.arange(len(resample_unique)), repeats=resample_count - 1)


        # -------
        this_win_start = this_win_end
        particles = [particles[idx] for idx in resample_unique] + [copy.deepcopy(particles[idx]) for idx in resample_repeated_idx_in_unique]

        """
        import pandas as pd
        df = pd.DataFrame ({'window_idx': this_win_idx,
                            'window_end': this_win_end,
                            'particle_no': range(params['n_SMC_particles']),
                            'n_infection': [this_particle.statevar[this_particle.last_idx, 3] for this_particle in particles],
                            'n_seg_before': [this_particle.n_segregating[this_win_idx, 0] for this_particle in particles],
                            'weight': this_w_vector[:, 0],
                            'after_sampling': k_matrix[this_win_idx, :] })

        this_win_start = this_win_end
        particles = [particles[idx] for idx in resample_unique] + [copy.deepcopy(particles[idx]) for idx in resample_repeated_idx_in_unique]        
        #--------
        df['n_seg_after'] = [this_particle.n_segregating[this_win_idx, 0] for this_particle in particles]
        with open(params['fname_particle_selection'], 'a') as f:
            df.to_csv(f, mode='a', header=f.tell() == 0, sep="\t")
        """
        #-------

    # calculate overall likelihood - if not "no_introduction"
    log_likelihood = 0
    for this_win_idx in range(len(sim_windows)):
        # sum if there is any data (for 1, 2, or 1 and 2)
        if (np.isnan(data['S_in_win'][this_win_idx])).all():    # no data exist for any of clades
            print (this_win_idx, ": no data")
            continue

        else:
            # calculate log likelihood
            w_matrix[this_win_idx, :] = np.where(w_matrix[this_win_idx, :] < 0, 0, w_matrix[this_win_idx, :])
            sum_likelihood_passage = sum(w_matrix[this_win_idx, :]) / params['n_SMC_particles']
            log_likelihood = log_likelihood + np.log(sum_likelihood_passage)

            printStr = '\t'.join([str(x) for x in [this_win_idx, data['data_flag'][this_win_idx, :], sum(w_matrix[this_win_idx, :]), sum_likelihood_passage, np.log(sum_likelihood_passage), log_likelihood]])
            """print(printStr, file=open(params['fname_likelihood_each_window'], 'a'))"""
            print(printStr)

    return log_likelihood, particles, False    # abandon t_intro




# ------------------
def _match_data_and_simwin(imported_data, sim_windows, n_clade):

    data = copy.deepcopy(imported_data)

    # todo - this could be simpler using numpy array
    # update data based on the simulation window; fill in missing win_interval, taj_d_in_win, var_taj_d_in_win, n_samples_in_win
    new_win_intervals = []
    new_taj_d_in_win = []
    new_var_taj_d_in_win = []
    new_n_samples_in_win = []

    data['win_intervals'] = np.array(data['win_intervals'])

    #print(max(data.win_intervals), params['t_start'], -1 * data.window_dt, "/", sim_windows)
    for win in sim_windows:
        if data['win_intervals'].tolist().count(win) == 0:  # if no corresponding data for this simulation window
            new_win_intervals.append(win)
            new_taj_d_in_win.append([np.nan] * n_clade)
            new_var_taj_d_in_win.append([np.nan] * n_clade)
            new_n_samples_in_win.append([0] * n_clade)

        else:  # if there is a data for this simulation window
            idx_in_data = data['win_intervals'].tolist().index(win)
            new_win_intervals.append(data['win_intervals'][idx_in_data])
            new_taj_d_in_win.append(data['S_in_win'][idx_in_data])
            new_n_samples_in_win.append(data['n_samples_in_win'][idx_in_data])

    data['win_intervals'] = np.array(new_win_intervals)
    data['S_in_win'] = np.array(new_taj_d_in_win)
    data['n_samples_in_win'] = np.array(new_n_samples_in_win)

    # new data flag type : after 2020-12-15
    data['data_flag'] = np.full_like(data['n_samples_in_win'], np.nan)
    before_first_data_point = -1
    between_first_and_last_data_point = 0
    after_last_data_point = 1

    for i in range(n_clade):

        win_intervals = data['win_intervals']
        n_samples_in_win = data['n_samples_in_win'][:, i]

        sample_point = np.multiply(n_samples_in_win > 0, ~(np.isnan(n_samples_in_win)))
        first_sample_point = win_intervals[sample_point].min()
        last_sample_point = win_intervals[sample_point].max()

        print("first_sample_point:", first_sample_point)
        print("last_sample_point:", last_sample_point)

        data_flag = np.where (win_intervals < first_sample_point, before_first_data_point, data['data_flag'][:, i])
        data_flag = np.where(np.multiply(first_sample_point <= win_intervals, win_intervals <= last_sample_point), between_first_and_last_data_point, data_flag)
        data_flag = np.where (last_sample_point < win_intervals, after_last_data_point, data_flag)

        # print (list(zip(n_samples_in_win, data_flag)))
        data['data_flag'][:, i] = data_flag


    return data


def _resample_particle (w, n_particle):
    # w = this_w_vector = weight of particles in this window
    # choose particles based on the weight
    if np.sum(w) == 0:
        # when all particles have gone extinct
        # or,  probTajD (in ‘run_SMC_SEIR’) gets evaluated to NaN for all particles,
            # which happens if the model predictions of all particles are exceptionally bad
            # or Tajima’s D can’t effectively be computed from the particles’ simulated data
        k = np.random.randint(low = n_particle, size = n_particle)

    else:
        w_norm = np.array(w) / np.sum(w)  # normalized weight -> now the sum is 1
        k = np.random.choice(range(n_particle), size=n_particle, replace=True, p=w_norm)

    return k
