# ------------------
## genotype_count.py (Salje_hetero; SE1E2IR model)
# ------------------
from numba import jit
import numpy as np
import random

import sys
sys.path.append('..')

import genotype_info as genotype_info_
from tools import misc_ as misc



# --------------
# model dependent functions
# 1. get the number of infection events
# 2. get the number of noninfection event & update
# 3. update the number of infection event
#    (so that the infection events and noninfection events do not interfere with each other)

def infection_event_mutate(count_event, statevar_curr, genotype_count, clade_idx, genotype_info, n_genotype, n_site, param_mu, rng):
    # model dependent

    src_idx = np.array([1, 2, 3, 4])
    dst_idx = np.array([0, 0, 0, 0])

    chosen_genotype_count_sum = np.zeros([len(dst_idx), n_genotype])
    n_new_genotype_sum = np.zeros([len(dst_idx)])

    n_genotype_unchanging = n_genotype
    for i in range(src_idx.size):
        s_idx = src_idx[i]

        n_ind_total = int(statevar_curr[s_idx+2])
        n_ind_to_pick = count_event[i]
        src_count = genotype_count[:n_genotype_unchanging, s_idx]

        if (n_ind_total != src_count.sum()):
            raise Exception(f'{s_idx}: n_ind_total {n_ind_total} is not equal to src_count.sum() {src_count.sum()}')

        if (n_genotype_unchanging > 0) and (n_ind_to_pick > 0):
            # choose genotype one infectious class and get the number of individuals with different mutation counts
            #print (src_count, n_ind_total, n_ind_to_pick, n_genotype_unchanging)
            chosen_genotype = rng.choice(n_genotype_unchanging, size=n_ind_to_pick, p=src_count/n_ind_total, replace=True)
            mutation_count = rng.multinomial(n=n_ind_to_pick, pvals=param_mu)

            chosen_genotype_count, n_new_genotype, genotype_info, n_site = _update_mutation(chosen_genotype, mutation_count, clade_idx, genotype_info, n_genotype, n_genotype_unchanging, n_site)

            chosen_genotype_count_sum[i] = chosen_genotype_count
            n_new_genotype_sum[i] = n_new_genotype

            n_genotype += n_new_genotype

            
    return chosen_genotype_count_sum, n_new_genotype_sum, genotype_info, n_site

@jit(nopython=True)
def infection_event_update(chosen_genotype_count, n_new_genotypes, p_genotype_count, clade_range, n_genotype):
    # model dependent

    src_idx = np.array([1, 2, 3, 4])
    dst_idx = np.array([0, 0, 0, 0])

    # print("chosen_genotype_count", chosen_genotype_count)
    if (n_genotype + n_new_genotypes.sum()) >= p_genotype_count.shape[0]:
        p_genotype_count = misc.add_row_to_nparray(p_genotype_count, n_genotype + n_new_genotypes.sum())

    n_genotype_orig = n_genotype
    for i in range(src_idx.size):
        s_idx = src_idx[i]
        d_idx = dst_idx[i]
        chosen_count = chosen_genotype_count[i]
        n_new_genotype = n_new_genotypes[i]

        # print(src_idx, ":", p_genotype_count[:, clade_range][:n_genotype, d_idx])
        p_genotype_count[:n_genotype_orig, d_idx] += chosen_count  # todo - multiclade
        n_genotype += n_new_genotype
        # print(src_idx, ":", p_genotype_count[:, clade_range][:n_genotype, d_idx])

        p_genotype_count[n_genotype_orig:n_genotype, d_idx] = 1  # todo - multiclade


    return p_genotype_count, n_genotype

def noninfection_event_update (count_event, statevar_curr, p_genotype_count, clade_range, n_genotype, rng):
    # model dependent

    src_idx = np.array([0, 1, 2, 0, 3, 4])[::-1]
    dst_idx = np.array([1, 2, 5, 3, 4, 5])[::-1]
    count_event = count_event[::-1]

    genotype_count = p_genotype_count[:, clade_range]

    for i in range(src_idx.size):
        s_idx = src_idx[i]
        d_idx = dst_idx[i]

        n_ind_total = int(statevar_curr[s_idx + 2])
        n_ind_to_pick = count_event[i]
        src_count = genotype_count[:n_genotype, s_idx]

        if (n_ind_total < src_count.sum()):     # when updating E-> I_h and then E->I_l, n_ind_total != src_count.sum() can happen
                                                # as statevar_curr is not updated yet but src_count is updated after E-> I_h
            raise Exception(
                f'{s_idx} , {d_idx}: n_ind_total {n_ind_total} is not equal to src_count.sum() {src_count.sum()}')

        if (src_count.sum() > 0) and (n_ind_to_pick > 0):
            chosen_ind = rng.choice(src_count.sum(), size=n_ind_to_pick, replace=False)
            chosen_genotype_count = _find_genotype_of_ind(chosen_ind, src_count, n_genotype).astype('int')
            _update_genotype_count(chosen_genotype_count, s_idx, d_idx, n_genotype, genotype_count)
            p_genotype_count[:, clade_range] = genotype_count

    return p_genotype_count


# --------------
# model independent functions

@jit(nopython=True)
def _update_mutation(chosen_genotype, mutation_count, clade_idx, genotype_info, n_genotype, n_genotype_unchanging, n_site):
    
    chosen_genotype_id = chosen_genotype[:mutation_count[0]]
    #print("chosen_genotype_id", chosen_genotype_id)
    #print ("np.bincount(genotype_list, minlength=n_genotype)", np.bincount(chosen_genotype_id, minlength=n_genotype_unchanging))

    chosen_genotype_count = _count_genotype(chosen_genotype_id, n_genotype_unchanging)
    #print("chosen_genotype_count", chosen_genotype_count)
    

    n_new_genotype = np.sum(mutation_count[1:])
    chosen_genotype_mu = chosen_genotype[mutation_count[0]:]
    genotype_info, n_site = genotype_info_.update_csr(clade_idx, chosen_genotype_mu, mutation_count[1:], genotype_info, n_genotype, n_site)

    return chosen_genotype_count, n_new_genotype, genotype_info, n_site

@jit(nopython=True)
def _count_genotype(genotype_list, n_genotype):
    return np.bincount(genotype_list, minlength=int(n_genotype))

@jit(nopython=True)
def _find_genotype_of_ind(chosen_ind, src_count, n_genotype):
    # genotype_count option 1 #
    ind_genotype = np.digitize(chosen_ind, bins=src_count.cumsum(), right=False);
    #print (">>ind_genotype:", ind_genotype, chosen_ind)
    chosen_genotype_count = np.bincount(ind_genotype, minlength=n_genotype)

    # genotype_count option 2 #chosen_genotype_count = np.histogram(chosen_ind, bins = src_count.cumsum(), range = (0, src_count.sum()))

    return chosen_genotype_count

@jit(nopython=True)
def _update_genotype_count(chosen_genotype_count, src_idx, dst_idx, n_genotype, genotype_count):
    genotype_count[:n_genotype, src_idx] -= chosen_genotype_count
    genotype_count[:n_genotype, dst_idx] += chosen_genotype_count

    return



