# calc_prior.py
## Author : Mike Martin


def calc_prior(params, theta):
    this_prior = 1
    # if there are no priors, will return 1
    if 'priors' in params.keys():
        priors = params['priors']
        for key, prior_func in priors.items():
            if type(key) is str:
                this_prior = \
                    this_prior * prior_func(theta[params['operators'].index(key)])
            elif type(key) is tuple:
                param_values = [theta[params['operators'].index(item)] for item in key]
                this_prior = \
                    this_prior * prior_func(*[theta[params.operators.index(item)] for item in key])
    return(this_prior)

