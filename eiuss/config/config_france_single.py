import numpy as np
from numba import jit
import importlib

# ------------------
def get_event_rate(params):
    ##################################################################
    # model specific values
    compartments = ['S', 'E1', 'E2_l', 'I_l', 'E2_h', 'I_h', 'R']
    infection_event = [['E2_l'  , 'E1'  , lambda x: x['beta'] * x['f_low']],
                       ['I_l'   , 'E1'  , lambda x: x['beta'] * x['f_low']],
                       ['E2_h'  , 'E1'  , lambda x: x['beta'] * x['f_high']],
                       ['I_h'   , 'E1'  , lambda x: x['beta'] * x['f_high']]]

    transition_event = [['E1'   , 'E2_l' , lambda x: x['gammaE1']],
                        ['E2_l' ,  'I_l' , lambda x: x['gammaE2']],
                        [ 'I_l' ,  'R'   , lambda x: x['gammaI']],
                        ['E1'   , 'E2_h' , lambda x: 0],  ## E1 -> E2_l and E2_h; binomial distribution
                        ['E2_h' ,  'I_h' , lambda x: x['gammaE2']],
                        [ 'I_h' ,  'R'   , lambda x: x['gammaI']]]
    ##################################################################

    ## infection events
    params['inv_infectious_period'] = 1/(1/params['gammaE2'] + 1/params['gammaI'])
    params['beta'] = (params['R0'] * params['inv_infectious_period']) / params['N']

    ## transmission heterogeneity
    params['c']      = ((1 - params['p_H']) / params['p_H']) / (1 / params['k'] - 1)
    params['f_low']  = (1 - params['k']) / (1 - params['p_H'])
    params['f_high'] = params['k'] / params['p_H']

    params['transition_split'] = [[0, 3, params['p_H']]]

    ## update parameters
    params['rate_infection']  = np.array([x[2](params) for x in infection_event])         ## x[2] = lambda function for each evemt rate
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
    param['rng2'] = param['rng']#np.random.default_rng(seed + 1)
    # --------------------
    #### DEFINE EPI MODEL ####
    ## model parameters ##
    param['mu']      = 0.33       # per birth mutation rate
    param['gammaE1'] = 1 / 4
    param['gammaE2'] = 1
    param['gammaI']  = 1 / 3
    param['R0']      = 2.5
    param['N']       = 67e6
    param['init']    = ['I_h']

    # transmission heterogeneity
    param['p_H']    = 0.10      # P(I_high) = r in Hao model
    param['k']      = 0.8       # proportion of infections caused by I_high

    # compartment model
    param['model_name']     = 'heteroSE1E2IR'
    param['compartments']   = [ 'E1', 'E2_l', 'I_l', 'E2_h', 'I_h', 'R']  ## except susceptible compartment
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

