import numpy as np
from scipy.stats import poisson


# ------------------
def get_weight_by_s(particles, window_data, first_datapoint, sample_class, n_grabs, rng):

    if window_data.size == 0:
        # if there is no data for this window, then weight depends only on extinct flag
        # this likely because we are earlier in the time series than where the data begins
        # or there is a missing data row and we're using homogeneous sampling
        # assumes we've appropriately trimmed empty data rows at the end of the time series
        print('no window data')
        w = []
        for particle_idx, particle in enumerate(particles):
            t_idx = particle.last_idx
            bool_E_sum_I = particle.genotype_count[:, -sample_class].sum() > 0

            if bool_E_sum_I == True:
                w.append(1)
            else:
                w.append(0)

        s = np.full(len(particles), np.nan)



    elif window_data['s'].isnull().values[0]:
        # if there is a data row but s but s is nan,
        # s is nan. implies this is a missing data point
        # thus, weight still depends only on extinction
        print('s is null')
        w = []
        for particle_idx, particle in enumerate(particles):
            t_idx = particle.last_idx
            bool_E_sum_I = particle.genotype_count[:, -sample_class].sum() > 0

            if bool_E_sum_I == True:
                w.append(1)
            else:
                w.append(0)

        s = np.full(len(particles), np.nan)



    elif window_data['n'].values[0] == 1 and window_data['s'].values[0] == 0:
        # if we are within the time series of the data,
        # but n = 1 and s = 0, then weight is determined by having at least 1
        # recovered individual during that windw
        print('only one sequence in time window')
        w = []
        for particle_idx, particle in enumerate(particles):
            t_idx = particle.last_idx
            bool_R = particle.genotype_count[:, sample_class].sum() > 0

            if bool_R == True:
                w.append(1)
            else:
                w.append(0)

        s = np.full(len(particles), np.nan)


    else:
        # finally, if there is actual data we can evaluate the likelihood
        # print(' > 1 sequence in time window')
        n_sample_data = int(window_data['n'].values[0])
        s = np.zeros(len(particles))

        for particle_idx, particle in enumerate(particles):
            segregating_site = []
            sampled_genotype_count = particle.genotype_count[:, sample_class]

            if sampled_genotype_count.sum() < n_sample_data:
                s[particle_idx] = np.nan

            else:
                for grab_idx in range(n_grabs):
                    chosen_sample_idx = rng.choice(sampled_genotype_count.sum(), size=n_sample_data, replace=False)
                    chosen_sample_genotype = np.digitize(chosen_sample_idx, bins=sampled_genotype_count.cumsum(), right=False)
                    segregating_site.append(
                        n_segregating_from_samples(chosen_sample_genotype.astype('int'), particle.genotype_info))

                s[particle_idx] = sum(segregating_site) / n_grabs

        w = poisson.pmf(window_data['s'], s)
        w = np.where(np.isnan(w), 0, w)
        print(f' > 1 sequence in time window, s = {s}, chosen_sample_genotype = {sampled_genotype_count.sum()}')

    # timestart values after earliest datapoint ahve a w of 0
    particle_t0 = np.array([particle.statevar[0, 0] for particle in particles])
    w = np.where(particle_t0 > first_datapoint, 0, w)


    return (w, s)






def n_segregating_from_samples (sample_genotypes, genotype_info):
    sample_genotypes = sample_genotypes.astype('int')

    col_idx = genotype_info[:, 0].astype('int')
    row_idx = genotype_info[:, 1].astype('int')

    sample_genotype_unique, sample_genotype_counts = np.unique(sample_genotypes, return_counts=True)

    n_sample_genotype_unique = len(sample_genotype_unique)
    n = sample_genotype_counts.sum()

    every_sites = np.array([], dtype='int')
    common_sites = col_idx[row_idx[sample_genotype_unique[0]]:row_idx[sample_genotype_unique[0]+1]]

    if n < 1:
        return np.nan

    else:
        if n_sample_genotype_unique == 1:
            return 0

        else:
            for i in range(n_sample_genotype_unique):
                genotype_i = sample_genotype_unique[i]
                sites_i = col_idx[row_idx[genotype_i]:row_idx[genotype_i + 1]]

                # segregating sites = all_sites - common_sites
                # ; a site is "segregating" unless it appears in all variants
                every_sites = np.union1d(every_sites, sites_i)
                common_sites = np.intersect1d(common_sites, sites_i)
            S = every_sites.size - common_sites.size
            return S





def run_smc (particles, params, imported_data, this_win_idx, this_win_start, this_win_end):

    no_introduction = False
    all_extinction = True

    first_datapoint = imported_data[imported_data['n'] >= 1]['window_end'].min()
    window_data     = imported_data[imported_data['window_end'] == this_win_end]

    for particle_idx in range(params['n_SMC_particles']):
        particle = particles[particle_idx]
        params['func_update_next_dt'](particle, params, this_win_start, this_win_end)       # update_next_dt

    this_w_vector, s = get_weight_by_s(particles, window_data, first_datapoint, params['idx_sampling'], params['n_grab'], params['rng'])         # default is 50 but we are trying 5

    for particle_idx, particle in enumerate(particles):
        particle.n_segregating[this_win_idx, 0] = s[particle_idx]



    return this_w_vector, no_introduction, all_extinction



