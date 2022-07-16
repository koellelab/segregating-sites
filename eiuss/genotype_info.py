# ------------------
## genotype_info.py 
## last modified : 2021-03-25 by YEONGSEON 
## just copied and pasted from [pmcmc_n_segregating_midas/_script/genotype_info.py]
# ------------------

# as genotype information changes only when new genotypes are added (mutation occurs),
# these functions will be used when picking new_S_to_E
# these functions are for each clade
import sys
import numpy as np
from numba import jit

import utils
# ------------------
@jit(nopython=True)
def update_csr_new (clade, mutation_parents, n_mutations, genotype_info, n_genotype, n_site):

    col_idx = genotype_info[:, clade*2 + 0]
    row_idx = genotype_info[:, clade*2 + 1]

    parent_unique = np.unique(mutation_parents)
    n_mutation_long = np.repeat(np.arange(1, n_mutations.size+1), n_mutations)

    if parent_unique.size == 0:
        return genotype_info, n_site


    for i in range(parent_unique.size):         # todo - fix this in midas version too
        parent = parent_unique[i]
        mutations = n_mutation_long[mutation_parents == parent]

        # find parent information
        parent_idx_start = row_idx[parent]  # the location where the parent column indicies starts
        parent_idx_end = row_idx[parent + 1]  # the location where the parent column indicies ends
        parent_sites = col_idx[parent_idx_start:parent_idx_end]
        n_site_parent = (parent_idx_end - parent_idx_start)  # the numnber of mutations(=non-zero data) in the parents


        # new row needed
        a = n_site_parent * (mutations.size) + mutations.sum()
        row_needed = (row_idx[n_genotype] + a) - genotype_info.shape[0] + 10
        if row_needed >= 0:
            new_arr = np.zeros((row_needed, genotype_info.shape[1]), dtype="int64")
            genotype_info = np.concatenate((genotype_info, new_arr), axis=0)

            col_idx = genotype_info[:, clade * 2 + 0]
            row_idx = genotype_info[:, clade * 2 + 1]


        # add new information
        for j in mutations:
            new_idx_start = row_idx[n_genotype]

            row_idx[n_genotype + 1] = row_idx[n_genotype] + n_site_parent + j
            col_idx[new_idx_start: new_idx_start + n_site_parent] = parent_sites
            col_idx[new_idx_start + n_site_parent: new_idx_start + n_site_parent + j] = np.arange(n_site, n_site + j)

            n_genotype += 1
            n_site += j

    return genotype_info, n_site





#@profile
@jit(nopython=True)
def update_csr (clade, mutation_parents, n_mutations, genotype_info, n_genotype, n_site):

    # Compressed sparse row (CSR, https://ko.wikipedia.org/wiki/%ED%9D%AC%EC%86%8C%ED%96%89%EB%A0%AC)
    # - row index: encoding the index in col_index where the given row starts (length of m+1)
    # - column index: column indicies of non-zero values (the length is the number of the nonzero entries in the original matrix)

    # array([[0, 0, 0, 0, 0, 0],
    #        [1, 0, 0, 0, 0, 0],
    #        [1, 1, 0, 0, 0, 0],
    #        [1, 0, 1, 1, 1, 0],
    #        [0, 0, 0, 0, 0, 0],
    #        [0, 0, 0, 0, 0, 0]])

    # => col_idx = [0 / 0, 1 / 0, 2, 3, 4]      # indic es
    #    row_idx = [0, 0, 1, 3, 7, 7, 7]        # indptr

    col_idx = genotype_info[:, clade*2 + 0]
    row_idx = genotype_info[:, clade*2 + 1]

    n_site_orig = n_site
    parent_idx = 0
    for idx, parent_counts in enumerate(n_mutations):
        i = idx + 1         # i = number of new mutation in an indiv / idx starts from 0
        parent_counts = n_mutations[idx]        # n_mutation starts from one mutation

        for j in range(parent_counts):
            parent = mutation_parents[parent_idx]
            parent_idx_start = row_idx[parent]  # the location where the parent column indicies starts
            parent_idx_end = row_idx[parent + 1]  # the location where the parent column indicies ends
            n_site_parent = (parent_idx_end - parent_idx_start)  # the numnber of mutations(=non-zero data) in the parents

            if (row_idx[n_genotype] + n_site_parent + i + 3 >= genotype_info.shape[0]):  # increase the array size if needed
                # print (genotype_info.shape, row_idx[n_genotype] + n_site_parent + i)
                genotype_info = utils.add_row_to_nparray(genotype_info, row_idx[n_genotype] + n_site_parent + i + 3)
                col_idx = genotype_info[:, clade * 2 + 0]
                row_idx = genotype_info[:, clade * 2 + 1]

            # 1. column index:
            # : this offspring will have n_site_parent (from parents) + i (new mutations) mutations
            new_idx_start = row_idx[n_genotype]         # "genotype n" is stored in (n-1)th row (as "genotype 0" is not stored)
            # -> information of "genotype n" is stored between row_idx[n_genotype-1]:row_idx[n_genotype]
            # -> information of "genotype n+1" is stored between row_idx[n_genotype]:row_idx[n_genotype+1]

            #print (parent, i, "|", new_idx_start, ",", n_genotype, parent, parent_idx_start, parent_idx_end, n_site_parent, "|",  new_idx_start +n_site_parent, "|", genotype_info.shape , row_idx[n_genotype] + n_site_parent + i)
            col_idx[new_idx_start: new_idx_start + n_site_parent] = col_idx[parent_idx_start:parent_idx_end]  # from parent
            col_idx[new_idx_start + n_site_parent: new_idx_start + n_site_parent + i] = np.arange(n_site, n_site + i)  # new mutations

            # 2. row index
            # : column index of this offspring ends at n_site_parent + i
            # : as there are n_site_parent + i non-zero data
            row_idx[n_genotype + 1] = row_idx[n_genotype] + n_site_parent + i

            n_genotype += 1
            n_site += i
            parent_idx += 1

    #if np.sum(n_mutations) != parent_idx:
    #    print (n_mutations, parent_idx)

    #if n_site_orig + np.sum(n_mutations * np.arange(1, n_mutations.size + 1)) != n_site:
    #if np.sum(n_mutations) > 0:
    #    print ("ya", n_site_orig, n_mutations, n_site)

    return genotype_info, n_site






@jit(nopython=True)
def update_csr_old (clade, mutation_parents, n_mutations, genotype_info ,n_genotype, n_site):

    # Compressed sparse row (CSR, https://ko.wikipedia.org/wiki/%ED%9D%AC%EC%86%8C%ED%96%89%EB%A0%AC)
    # - row index: encoding the index in col_index where the given row starts (length of m+1)
    # - column index: column indicies of non-zero values (the length is the number of the nonzero entries in the original matrix)

    # array([[0, 0, 0, 0, 0, 0],
    #        [1, 0, 0, 0, 0, 0],
    #        [1, 1, 0, 0, 0, 0],
    #        [1, 0, 1, 1, 1, 0],
    #        [0, 0, 0, 0, 0, 0],
    #        [0, 0, 0, 0, 0, 0]])

    # => col_idx = [0 / 0, 1 / 0, 2, 3, 4]
    #    row_idx = [0, 0, 1, 3, 7, 7, 7]

    col_idx = genotype_info[:, clade*2 + 0]
    row_idx = genotype_info[:, clade*2 + 1]


    for idx, i in enumerate(n_mutations):
        parent = mutation_parents[idx]
        #print ("hey1")
        parent_idx_start = row_idx[parent]  # the location where the parent column indicies starts
        parent_idx_end = row_idx[parent+1]        # the location where the parent column indicies ends
        n_site_parent = (parent_idx_end - parent_idx_start)     # the numnber of mutations(=non-zero data) in the parents

        #print(">>>", i, genotype_info.shape, n_genotype, parent,  "|",parent_idx_start, parent_idx_end,  "|",row_idx[n_genotype] , n_site_parent, i, row_idx[n_genotype] + n_site_parent + i + 1)

        if ( row_idx[n_genotype] + n_site_parent + i + 1 >= genotype_info.shape[0]):         # increase the array size if needed
            #print (genotype_info.shape, row_idx[n_genotype] + n_site_parent + i)
            genotype_info = utils.add_row_to_nparray(genotype_info, row_idx[n_genotype] + n_site_parent + i + 1)
            col_idx = genotype_info[:, clade * 2 + 0]
            row_idx = genotype_info[:, clade * 2 + 1]
            #print ("genotype_array_added")

        #print("hey2")
        # 1. column index:
        # : this offspring will have n_site_parent (from parents) + i (new mutations) mutations
        new_idx_start = row_idx[n_genotype]         # "genotype n" is stored in (n-1)th row (as "genotype 0" is not stored)
        # -> information of "genotype n" is stored between row_idx[n_genotype-1]:row_idx[n_genotype]
        # -> information of "genotype n+1" is stored between row_idx[n_genotype]:row_idx[n_genotype+1]

        #print (parent, i, "|", new_idx_start, ",", n_genotype, parent, parent_idx_start, parent_idx_end, n_site_parent, "|",  new_idx_start +n_site_parent, "|", genotype_info.shape , row_idx[n_genotype] + n_site_parent + i)
        col_idx[new_idx_start : new_idx_start +n_site_parent] = col_idx[parent_idx_start:parent_idx_end]                         # from parent
        #print("hey4-1")
        col_idx[new_idx_start + n_site_parent : new_idx_start + n_site_parent + i] = np.arange(n_site, n_site + i)    # new mutations
        #print("hey4-2")

        # 2. row index
        # : column index of this offspring ends at n_site_parent + i
        # : as there are n_site_parent + i non-zero data
        row_idx[n_genotype+1] = row_idx[n_genotype] + n_site_parent + i

        n_genotype += 1
        n_site += i
        #print("hey5")

    #print ("><", n_genotype)

    return genotype_info, n_genotype, n_site
    #return genotype_info, n_site_
