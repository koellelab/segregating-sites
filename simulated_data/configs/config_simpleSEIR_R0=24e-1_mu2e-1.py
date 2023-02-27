import numpy as np
from numba import jit
import importlib

# ------------------
def get_event_rate(params):
    infection_event   = [['I', 'E', lambda x: x['beta']]]
    transition_event  = [['E', 'I', lambda x: x['gammaE']],
                         ['I', 'R', lambda x: x['gammaI']]]

    ## infection events
    params['inv_infectious_period'] = params['gammaI']
    params['beta'] = (params['R0'] * params['inv_infectious_period']) / params['N']
    params['rate_infection'] = np.array([x[2](params) for x in infection_event])        ## x[2] = lambda function for each evemt rate

    ## transmision_event
    params['rate_transition'] = np.array([x[2](params) for x in transition_event])        ## x[2] = lambda function for each evemt rate

    return params


def update_dependent_param(param):
    from scipy.stats import poisson
    param['mu_multinom'] = np.array([poisson.pmf(x, param['mu']) for x in range(6)])
    param = get_event_rate(param)

    return param




def set_param(seed):
    param = {}
    param['seed'] = seed
    param['rng']  = np.random.default_rng(seed)
    param['rng2'] = np.random.default_rng(seed + 1)
    # --------------------
    #### DEFINE EPI MODEL ####
    ## model parameters ##
    param['mu']      = 0.2       # per birth mutation rate
    param['R0']      = 2.4
    param['gammaE']  = 1 / 2
    param['gammaI']  = 1 / 3
    param['N']       = 1e5
    param['init']    = ['I']


    #param['n_clades'] = 1
    # compartment model
    param['model_name']     = 'simpleSEIR'
    param['compartments']   = ['E', 'I', 'R']    ## except susceptible compartment
    param['n_class']        = len(param['compartments'])

    compartments = {x: idx for idx, x in enumerate(param['compartments'])}
    param['idx_sampling']   = compartments['R']
    param['init']           = np.isin(param['compartments'], param['init']).astype('int')

    model_nextdt            = importlib.import_module("models." + param["model_name"] + ".next_dt")
    param['model_nextdt']   = model_nextdt.update
    # --------------------
    #### TIME SERIES CONFIGURATION ####
    param['timestart']      = 0
    param['timeend']        = 200
    param['dt']             = 0.1  # time window for gillespie simulations

    # ====================
    #### BELOW THIS DO NOT EDIT ####
    # calculate some parameters that depend on other parameters
    param = update_dependent_param(param)
    print ("------------------------------")
    for key, value in param.items():
        print (f'{key:>15s} : {value}')#            (key, "\t:", value)
    print("------------------------------")

    return param
