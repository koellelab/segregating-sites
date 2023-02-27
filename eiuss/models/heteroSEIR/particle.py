import numpy as np
# ------------------
class Particle():
    def __init__(self, n_clade, n_class, init_size = 10):
        ## numpy array allocation
        self.n_clade = 1
        self.n_class = n_class

        self.statevar       = np.zeros([init_size, n_class + 2], dtype='float64')    ## t | S | [other compartments]
        self.genotype_count = np.zeros([init_size, n_class], dtype='int64')        ## [other compartments]
        self.genotype_info  = np.zeros([init_size, 2], dtype='int64')

        self.n_segregating = np.full([init_size, 2], np.nan, dtype='float64')  ## segregating site | window_end

        self.n_genotype     = np.array([0])
        self.n_site         = np.array([0])
        self.last_idx       = 0





    def initiate (self, param, win_length):
        ## set initial values for state variables
        self.statevar[0, 0]       = np.min(param['timestart'])
        self.statevar[0, 1]       = param['N'] - param['init'].sum()
        self.statevar[0, 2:]      = param['init']
        self.genotype_count[0, :] = param['init']


        self.n_genotype[0] = 1
        self.n_site    [0] = 0
        self.last_idx      = 0

        return



    def trim_arrays(self):
        self.statevar       = self.statevar[:self.last_idx+1, :]
        self.genotype_count = self.genotype_count[:np.max(self.n_genotype), :]

        return