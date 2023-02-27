# ------------------
## statistics.py 
## last modified : before 2021-03-25 by YEONGSEON 
## just copied and pasted from [pmcmc_n_segregating_midas/_script/weight_by_S/statistics.py]
# ------------------
import numpy as np


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


def n_segregating_from_samples_faster (sample_genotypes, genotype_info):
    sample_genotypes = sample_genotypes.astype('int')

    col_idx = genotype_info[:, 0].astype('int')
    row_idx = genotype_info[:, 1].astype('int')

    sample_genotype_unique, sample_genotype_counts = np.unique(sample_genotypes, return_counts=True)
    n_sample_genotype_unique = len(sample_genotype_unique)
    n = sample_genotype_counts.sum()

    if n < 1:
        return np.nan
    else:
        if n_sample_genotype_unique == 1:
            return 0
        else:
            every_sites = []
            for i in range(n_sample_genotype_unique):
                genotype_i = sample_genotype_unique[i]
                sites_i = col_idx[row_idx[genotype_i]:row_idx[genotype_i + 1]]

                # segregating sites = all_sites - common_sites
                # ; a site is "segregating" unless it appears in all variants
                this_genotype_count = sample_genotype_counts[i]
                every_sites += sites_i.tolist() * this_genotype_count

            site_frequency = np.unique(every_sites, return_counts=True)
            S = (site_frequency[1] < n).sum()

            return S
