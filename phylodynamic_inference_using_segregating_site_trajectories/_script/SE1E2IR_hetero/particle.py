# ------------------
## particle.py (Salje_hetero; SE1E2IR model)
# ------------------
import numpy as np
#sys.path.append('..')

from tools import misc_ as misc
# ------------------
class Particle():
    def __init__(self, n_clade, n_class, init_size = 10):
        self.n_clade = n_clade  # either from data or given
        self.n_class = n_class        # E1, E2, I, cumI_or_R (except S class)

        self.statevar = np.empty([init_size, 2 + self.n_class * n_clade], dtype='float64')     # t should be in float type
        self.genotype_count = np.empty([init_size, self.n_class * n_clade], dtype='int64')
        self.genotype_info = np.empty([init_size, 2 * n_clade], dtype = 'int64')            # there will be only row_idx and col_idx for each clade

        self.n_segregating = np.empty([init_size, n_clade+1], dtype='float64')  # the last column is t_for_tajD

        self.n_genotype = np.empty(n_clade, dtype='int64')
        self.n_site = np.empty(n_clade, dtype='int64')
        self.last_idx = 0



    def introduce_clade(self, t_idx, clade_introduction, idx_intro_s, idx_intro_g):
        # model dependent
        clade_to_update = np.repeat(clade_introduction, self.n_class)

        ## 1. state_var - high introduced first
        state_var_to_update = np.tile([True, False, False, False, False, True],self.n_clade)  # E introduced first
        self.statevar[t_idx, 2:][(state_var_to_update & clade_to_update)] = 1                 # increase introduced class by 1
        self.statevar[t_idx, 1] = self.statevar[t_idx, 1] - 1                                 # reduce S by 1

        ## 2. genotype_counts - high introduced first
        genotype_count_to_update = np.tile([True, False, False, False, False, False],self.n_clade)  # E introduced first
        self.n_genotype[clade_introduction] = 1
        self.genotype_count[0, (genotype_count_to_update & clade_to_update)] = 1

        return


    def initiate(self, params, sim_win_len):
        if (sim_win_len >= self.n_segregating.shape[0]):
            self.n_segregating = misc.add_row_to_nparray(self.n_segregating, sim_win_len + 1)

        self.statevar[:] = 0
        self.genotype_count[:] = 0
        self.genotype_info[:] = 0

        self.n_segregating[:] = np.nan

        self.n_genotype[:] = 0
        self.n_site[:] = 0
        self.last_idx = 0


        # introduce the very first clade
        if 'timeintro' in params.keys():
            starting_clade = (params['timeintro'] == np.min(params['timeintro']))
        else :
            starting_clade = [True]
        self.statevar[0, 0] = np.min(params['timeintro'])
        self.statevar[0, 1] = params['N']

        self.introduce_clade(0, starting_clade, params['idx_clade_intro_statevar'], params['idx_clade_intro_genocount'])  # reduce S_curr by 1 and set I_init and cumI_init of the starting clade

        return


    def trim_arrays(self):
        self.statevar = self.statevar[:self.last_idx+1, :]
        self.genotype_count = self.genotype_count[:np.max(self.n_genotype), :]

        return


