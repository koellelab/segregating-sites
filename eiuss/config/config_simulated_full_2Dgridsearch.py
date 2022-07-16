import os
import numpy as np
from numba import jit
import importlib


import utils
# ------------------
@jit(nopython=True)
def event_rate_infection(param_R0, param_N, param_gammaI, n_clade):
    return np.array([param_R0 * param_gammaI])/ param_N


@jit(nopython=True)
def event_rate_noninfection(param_gammaE, param_gammaI, n_clade):
    return np.array([param_gammaE, param_gammaI])   # event_rate for high group is 0 -> see assign_high_group function




def get_event_rate(params):
    params['rate_noninfection'] = event_rate_noninfection(params['gammaE'], params['gammaI'],params['n_clades'])
    params['rate_infection'] = event_rate_infection(params['R0'], params['N'], params['gammaI'], params['n_clades'])
    return params




def set_param(imported_data):
    param = {}
    param['rng'] = np.random.default_rng()  # random number generator
    # --------------------
    #### DEFINE EPI MODEL ####
    ## model parameters ##
    # composite parametes must be defined as function
    # composite parameters must be defined in order
    param['mu'] = 0.2           # per birth mutation rate
    param['gammaE'] = 1 / 2       # calculated with supplementary_info/distribution_Salje_E_class/salje_distribution.py
    param['gammaI'] = 1 / 3
    param['R0'] = 1.6

    # initial conditions
    param['N'] = 1e5
    param['n_clades'] = 1

    # compartment model
    compartments            = ['E', 'I', 'R']
    param['n_class']        = len(compartments)    
    param['idx_sampling']   = 2
    
    param['model_name']          = 'model_SEIR'
    model_nextdt                 = importlib.import_module(param["model_name"] + ".next_dt")
    param['func_update_next_dt'] = model_nextdt.update

    # --------------------
    #### TIME SERIES CONFIGURATION ####
    param['timestart']      = 0
    param['timeend']        = max(imported_data['window_end'])  # date in matlab ordinal format for simulations to end
    param['win_intervals']  = imported_data['window_end']  # dates of window intervals
    param['window_dt']      = imported_data['window_end'][1] - imported_data['window_end'][0]
    param['dt']             = 0.1  # time window for gillespie simulations
    # --------------------
    #### SMC CONFIGURATION ####
    param['n_SMC_particles'] = 200  # number of particles
    param['n_grab']          = 50
    #### INFERENCE CONFIGURATION ####
    param['inference_framework'] = 'gridsearch.py'  # python script which holds the inference method (gridsearch, mcmc, mif, etc.)
    param['n_reps']              = 20
    param['operators']           = ['R0', 'timestart']  # must be key in param or param['model_params'] dictionaries
    param['range']               = [utils.np_arrange_with_endpoint(1.00, 2.50, 0.1), np.arange(-50, 40, 2)]
    
    param['out_name']            = "/simulated_data/prop500_full/simulated_data_prop500_full_2Dgridsearch"
    # --------------------
    print("------------------------------")
    for key, value in param.items():
        print(key, "\t:", value)
    print("------------------------------")


    #### BELOW THIS DO NOT EDIT ####
    # auto calculated based on other parameters
    # precalcualte some rates
    from scipy.stats import poisson
    param['mu_multinom'] = np.array([poisson.pmf(x, param['mu']) for x in range(6)])  # 0 mutation ~ 5 mutations
    param = get_event_rate(param)

    return param


## grep 'profilelikelihood_1028' config_simulated_2Dprofilelikelihood*
## grep ' # range of values to explore' config_simulated_2Dprofilelikelihood*
## grep 'tsv' config_simulated_2Dprofilelikelihood*
##sed -i'.original' -e 's|[mean_logL] piuss_simulated1234-500samples_gridsearch_1028|[mean_logL] piuss_simulated1234-500samples_gridsearch_1028|g' conf*2Dprofilelikelihood*.py 




