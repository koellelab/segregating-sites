import numpy as np
from scipy.stats import poisson

def check_and_fill_in_params (params_pmcmc):

    # [check errors] time-related parameters that depend on other time-related parameters
    if np.any(params_pmcmc['timestart'] > params_pmcmc['timeintro']):
        raise Exception("ParameterError: time_start should be earlier (smaller) than any time_intro")
    if np.any(params_pmcmc['timestart'] > params_pmcmc['timeend']):
        raise Exception("ParameterError: time_end should be later (largerer) than time_start")
    if np.any(params_pmcmc['timeintro'] > params_pmcmc['timeend']):
        raise Exception("ParameterError: time_end should be later (largerer) than any param_timeintro")


    # [check errors] operators (which parameters are to be operated on in MCMC chain)
    for item in params_pmcmc['operators']:
        if item not in params_pmcmc.keys():
            raise Exception(f'cannot operate on parameter {item} because it is not defined')

    # [check errors] SD of the proposal distribution for each operator
    for item in params_pmcmc['operator_SD'].keys():
        if item not in params_pmcmc.keys():
            raise Exception(f'cannot set operator SD on parameter {item} because it is not defined')

    # [check errors] test if all keys are part of the operators list
    for key in params_pmcmc['priors']:
        if type(key) is str:
            if key not in params_pmcmc['operators']:
                raise Exception(f'item {key} in priors dictionary not in operators list')
        elif type(key) is list or type(key) is tuple:
            for item in key:
                if key not in params_pmcmc['operators']:
                    raise Exception(f'item {item} in priors dictionary not in operators list')
        else:
            raise Exception(f'keys in priors dictionary must be of type: string, list, or tuple')


    # ------------------
    params_pmcmc['theta_init'] = [params_pmcmc[key] for key in params_pmcmc['operators']]
    params_pmcmc['theta_sd'] = [params_pmcmc['operator_SD'][key] for key in params_pmcmc['operators']]
    params_pmcmc['theta_sd'] = np.diag(params_pmcmc['theta_sd'])
    params_pmcmc['theta_var'] = params_pmcmc['theta_sd'] ** 2



    return params_pmcmc




def check_time_params (params_pmcmc):
    # [check errors] time-related parameters that depend on other time-related parameters
    if np.any(params_pmcmc['timestart'] > params_pmcmc['timeintro']):
        raise Exception("ParameterError: time_start should be earlier (smaller) than any time_intro")
    if np.any(params_pmcmc['timestart'] > params_pmcmc['timeend']):
        raise Exception("ParameterError: time_end should be later (largerer) than time_start")
    if np.any(params_pmcmc['timeintro'] > params_pmcmc['timeend']):
        raise Exception("ParameterError: time_end should be later (largerer) than any param_timeintro")

    return params_pmcmc


def check_pmcmc_params (params_pmcmc):
    # [check errors] operators (which parameters are to be operated on in MCMC chain)
    for item in params_pmcmc['operators']:
        if item not in params_pmcmc.keys():
            raise Exception(f'cannot operate on parameter {item} because it is not defined')

    # [check errors] SD of the proposal distribution for each operator
    for item in params_pmcmc['operator_SD'].keys():
        if item not in params_pmcmc.keys():
            raise Exception(f'cannot set operator SD on parameter {item} because it is not defined')

    # [check errors] test if all keys are part of the operators list
    for key in params_pmcmc['priors']:
        if type(key) is str:
            if key not in params_pmcmc['operators']:
                raise Exception(f'item {key} in priors dictionary not in operators list')
        elif type(key) is list or type(key) is tuple:
            for item in key:
                if key not in params_pmcmc['operators']:
                    raise Exception(f'item {item} in priors dictionary not in operators list')
        else:
            raise Exception(f'keys in priors dictionary must be of type: string, list, or tuple')


    # ------------------
    params_pmcmc['theta_init'] = [params_pmcmc[key] for key in params_pmcmc['operators']]
    params_pmcmc['theta_sd'] = [params_pmcmc['operator_SD'][key] for key in params_pmcmc['operators']]
    params_pmcmc['theta_sd'] = np.diag(params_pmcmc['theta_sd'])
    params_pmcmc['theta_var'] = params_pmcmc['theta_sd'] ** 2



    return params_pmcmc