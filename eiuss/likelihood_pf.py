import copy
import numpy as np
import pandas as pd

import smc
# ------------------
def get_statevars(particles, i, win_idx):
    t_idx_win = particles[i].last_idx
    state_win = particles[i].statevar[t_idx_win, :]
    seg_win = particles[i].n_segregating[win_idx, 0]

    return state_win.tolist() + [seg_win]


def calculate (particles, params, imported_data, config):
    ## update dependent parameters
    params = config.update_dependent_param(params)

    # get simulation window
    first_datapoint = imported_data[imported_data['n'] >= 1]['window_end'].min()
    if first_datapoint < params['timestart']:
        return -np.inf, None, None

    sim_windows = np.arange(max(imported_data['window_end']), params['timestart'], params['window_dt'] * (-1))
    sim_windows = sim_windows[::-1]

    print("------------------------------")
    print(">> updated parameters and dependent parameters")
    for key, value in params.items():
        if key in params['operators'] + ['rate_infection', 'rate_transition','mu_multinom']:
            print(f'{key:>15s} : {value}')  # (key, "\t:", value)

        if key == "eta":
            print(f'{key:>15s} : {value}')  # (key, "\t:", value)
    print("------------------------------")


    # initialize particles & add n_segregating row if needed
    for particle in particles:
        particle.initiate(params, len(sim_windows))
        particle.n_segregating[:len(sim_windows), -1] = sim_windows

    # matrix for weight, resampling
    w_matrix = np.full([len(sim_windows), params['n_SMC_particles']], np.nan)                 # shape = sim_window * n_SMC_particles
    k_matrix = np.full([len(sim_windows), params['n_SMC_particles']], -9999, 'int')           # shape  = sim_window * n_SMC_particles


    col_names = "t_before s_before e_before i_before r_before seg_before t_after s_after e_after i_after r_after seg_after selected".split(" ")
    df_particles = pd.DataFrame()


    # simulate and resample particles for each window
    this_win_start = params['timestart']
    for this_win_idx in range(len(sim_windows)):

        saved_particle = [get_statevars(particles, i, this_win_idx-1) for i in range(params['n_SMC_particles'])]       ## log how particles evolved

        this_win_end = sim_windows[this_win_idx]
        print("loglikelihood :", (this_win_idx, this_win_start, this_win_end), flush=True)
        this_w_vector = smc.run_smc(particles, params, imported_data, this_win_idx, this_win_start, this_win_end)

        w_matrix[this_win_idx, :] = this_w_vector
        k_matrix[this_win_idx, :] = _resample_particle(np.where(this_w_vector < 0, 0, this_w_vector), params['n_SMC_particles'])

        ## log how particles evolved
        saved_particle = [x+get_statevars(particles, i, this_win_idx) for i, x in enumerate(saved_particle)]
        saved_particle = [x+["V"] if np.isin(i, k_matrix[this_win_idx, :]) else x+["X"] for i, x in enumerate(saved_particle)]


        # resample particles based on their weight
        resample_unique, resample_count = np.unique(k_matrix[this_win_idx, :], return_counts=True)
        resample_repeated_idx_in_unique = np.repeat(resample_unique, repeats=resample_count - 1)
        particles = [particles[idx] for idx in resample_unique] + [copy.deepcopy(particles[idx]) for idx in resample_repeated_idx_in_unique]

        this_win_start = this_win_end

        df_particles = pd.concat([df_particles,
                                  pd.DataFrame(np.array(saved_particle), columns=col_names)])

        #print(df_particles.to_string())


    # calculate overall likelihood
    data_exists = \
        np.where(np.isin(sim_windows, imported_data[~imported_data['s'].isnull()]['window_end']) == True)[0]

    # replace nans with 0s
    w_matrix = np.where(~(w_matrix >= 0), 0, w_matrix)

    # average likelihoods over particles
    lh = w_matrix[data_exists, :].sum(axis=1) / params['n_SMC_particles']

    # log likelihood score for this parameter set
    llh = np.log(lh)
    llh_sum = llh.sum()

    print(pd.DataFrame(np.vstack([data_exists, sim_windows[data_exists], lh,llh]).T,
                       columns=['window', 'window_end','lh', 'llh']).to_string())


    df_particles.to_csv(params['out_name'] + '.tsv', sep="\t", index=False)


    return llh_sum, particles, lh


# ------------------
def _resample_particle (w, n_particle):
    # w = this_w_vector = weight of particles in this window
    # choose particles based on the weight
    if np.sum(w) == 0:
        # when all particles have gone extinct
        # or, weight  gets evaluated to NaN for all particles,
            # which happens if the model predictions of all particles are exceptionally bad
            # or the number of segregating sites  can’t effectively be computed from the particles’ simulated data
        k = np.random.randint(low = n_particle, size = n_particle)

    else:
        w_norm = np.array(w) / np.sum(w)  # normalized weight -> now the sum is 1
        k = np.random.choice(range(n_particle), size=n_particle, replace=True, p=w_norm)

    return k
