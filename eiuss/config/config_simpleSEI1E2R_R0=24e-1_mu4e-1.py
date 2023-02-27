import numpy as np
from numba import jit




def compartment_model (param):
    param['model_name']   = 'hSEIR'
    param['compartments'] = ['E', 'I_l', 'I_h', 'R']
    param['n_class']      = len(param['compartments'])

    infectious_period =  1 / param['gammaI']
    infection_event   = [['I_l', 'E', (param['R0'] * (1/infectious_period) * param['f_low'])  / param['N']],
                         ['I_h', 'E', (param['R0'] * (1/infectious_period) * param['f_high']) / param['N']]]

    transition_event = [['E'  , 'I_l', param['gammaE']],          ## E1 -> I_l, I_h
                        ['E'  , 'I_h', 0              ],
                        ['I_l', 'R'  , param['gammaI']],
                        ['I_h', 'R'  , param['gammaI']]]
    param['transition_split'] = [[[0, 1], (1 - param['p_H'], param['p_H'])]]

    sampling_compartment = ['R']
    # -----------------
    compartments = {x: idx for idx, x in enumerate(param['compartments'])}
    param['infection_S']    = np.array([-1])
    param['infection_cumI'] = np.array([param['n_class']])
    param['infection_src' ] = np.array([compartments[x[0]] for x in infection_event])
    param['infection_dst' ] = np.array([compartments[x[1]] for x in infection_event])
    param['infection_rate'] = np.array([x[2]               for x in infection_event])

    param['transition_src' ] = np.array([compartments[x[0]] for x in transition_event])
    param['transition_dst' ] = np.array([compartments[x[1]] for x in transition_event])
    param['transition_rate'] = np.array([x[2]               for x in transition_event])

    param['sampling_src']    = np.array([compartments[x] for x in sampling_compartment])

    param['init'] = [param['compartment_init'][x] if x in param['compartment_init'].keys() else 0 for x in param['compartments']]
    param['init'] = np.array([param['N'] - sum(param['init'])] + param['init'])

    # -----------------
    param['rate_noninfection'] = param['transition_rate']
    param['rate_infection']    = param['infection_rate']
    param['idx_g_R']           = param['sampling_src']

    E = 0;
    I = 1;
    R = 2
    n_class = param['n_class']
    param['n_clades'] = 1
    param['idx_clade_intro_statevar'] = np.isin(np.arange(n_class) % n_class, [I, R])      # highly transmissible infectious individual is introduced
    param['idx_clade_intro_genocount'] = np.isin(np.arange(n_class) % n_class, [I])


    return param






def set_param_simul(seed):
    param = {}
    param['seed'] = seed
    param['rng']  = np.random.default_rng(seed)
    param['rng2'] = param['rng']
    # --------------------
    #### DEFINE EPI MODEL ####
    ## model parameters ##
    param['mu']      = 0.4       # per birth mutation rate
    param['R0']      = 2.4
    param['gammaE']  = 1 / 2
    param['gammaI']  = 1 / 3
    param['N']       = 1e5
    param['compartment_init'] = {'I_h': 1}

    # transmission heterogeneity
    param['p_H'] = 0.80  # P(I_high) = r in Hao model
    param['k'] = 0.8  # proportion of infections caused by I_high
    param['c'] = ((1 - param['p_H']) / param['p_H']) / (1 / param['k'] - 1)
    param['f_low'] = (1 - param['k']) / (1 - param['p_H'])
    param['f_high'] = param['k'] / param['p_H']

    # --------------------
    #### TIME SERIES CONFIGURATION ####
    param['timestart']      = 0
    param['timeend']        = 200
    param['dt']             = 0.1  # time window for gillespie simulations
    param['window_dt']      = 4
    # --------------------
    param['timeintro'] = 0

    # ====================
    #### BELOW THIS DO NOT EDIT ####
    param = compartment_model(param)
    # calculate some parameters that depend on other parameters
    from scipy.stats import poisson
    param['mu_multinom'] = np.array([poisson.pmf(x, param['mu']) for x in range(6)])




    print ("------------------------------")
    for key, value in param.items():
        print (key, "\t:", value)
    print("------------------------------")

    return param