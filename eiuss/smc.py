import numpy as np
import pandas as pd
from scipy.stats import poisson

import statistics_S
# -----------------
def get_weight_by_s(particles, window_data, sample_class, n_grabs, rng):

    s = np.full(len(particles), np.inf)
    w = np.full(len(particles), np.nan)


    if window_data.size == 0:
        # if there is no data for this window, then weight depends only on extinct flag
        # this likely because we are earlier in the time series than where the data begins
        # or there is a missing data row and we're using homogeneous sampling
        # assumes we've appropriately trimmed empty data rows at the end of the time series
        print('no window data')
        for particle_idx, particle in enumerate(particles):
            last_statevar = particle.statevar[particle.last_idx, :]
            s[particle_idx] = np.nan
            w[particle_idx] = (last_statevar[2:sample_class+2].sum() > 0).astype('int')       ## check not extinction
            #print("statevar:",last_statevar.tolist(),  last_statevar[2:sample_class + 2].tolist())

    elif window_data['n'].values[0] == 0:
        for particle_idx, particle in enumerate(particles):
            last_statevar = particle.statevar[particle.last_idx, :]
            s[particle_idx] = np.nan
            w[particle_idx] = (last_statevar[2:sample_class+2].sum() > 0).astype('int')       ## check not extinction
            #print("statevar:",last_statevar.tolist(),  last_statevar[2:sample_class + 2].tolist())

    elif window_data['n'].values[0] == 1:
        # if we are within the time series of the data,
        # but n = 1 and s = 0, then weight is determined by having at least 1
        # recovered individual during that windw
        if window_data['s'].values[0] != 0:
            raise Exception("Somethings wrong with the input data: S couldn't be other than 0 if n = 1")

        for particle_idx, particle in enumerate(particles):
            last_statevar = particle.statevar[particle.last_idx, :]
            s[particle_idx] = 0 if (particle.genotype_count[:, sample_class].sum() > 0) else np.nan
            w[particle_idx] = (particle.genotype_count[:, sample_class].sum() > 0).astype('int')       ## check_there_is_recovered

            #print(particle.genotype_count[:, :].sum(axis=0).tolist(), particle.genotype_count[:, sample_class].sum() )

    else:
        # finally, if there is actual data we can evaluate the likelihood
        # print(' > 1 sequence in time window')
        from scipy.stats import poisson
        n_sample_data = int(window_data['n'].values[0])

        for particle_idx, particle in enumerate(particles):
            #print(particle.genotype_count[:, :].sum(axis=0).tolist(), particle.genotype_count[:, sample_class].sum() )
            segregating_site = []
            sampled_genotype_count = particle.genotype_count[:, sample_class]

            if sampled_genotype_count.sum() < n_sample_data:
                s[particle_idx] = np.nan
                w[particle_idx] = 0
            else:
                for grab_idx in range(n_grabs):
                    chosen_sample_idx = rng.choice(sampled_genotype_count.sum(), size=n_sample_data, replace=False)
                    chosen_sample_genotype = np.digitize(chosen_sample_idx, bins=sampled_genotype_count.cumsum(), right=False)

                    segregating_site.append(
                        statistics_S.n_segregating_from_samples_faster(chosen_sample_genotype.astype('int'), particle.genotype_info))

                s[particle_idx] = sum(segregating_site) / n_grabs
                w[particle_idx] = poisson.pmf(window_data['s'], s[particle_idx])

        if np.isnan(w).any():
            print(s[np.where(np.isnan(w))], w[np.where(np.isnan(w))])
            raise Exception ("w = np.nan")
        w = np.where(np.isnan(w), 0, w)

    return (w, s)


def run_smc (particles, params, imported_data, this_win_idx, this_win_start, this_win_end):

    for idx in range(params['n_SMC_particles']):
        particle = particles[idx]
        params['model_nextdt'](particle, params, this_win_start, this_win_end)       # update_next_dt

    ## get weight by comparing with data
    window_data      = imported_data[imported_data['window_end'] == this_win_end]
    this_w_vector, s = \
        get_weight_by_s(particles, window_data, params['idx_sampling'], params['n_grabs'], params['rng'])

    ## update n_segregating for each particle
    for particle_idx, particle in enumerate(particles):
        particle.n_segregating[this_win_idx, 0] = s[particle_idx]

    return this_w_vector






















