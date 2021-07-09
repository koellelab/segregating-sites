# ------------------
## weight.py 
# ------------------
import numpy as np
from scipy.stats import poisson
from scipy.special import factorial


def get_weight_by_S (bool_cumI, bool_E_sum_I, data_flag, particle_S, data_S, data_tajD_var):

    # bool_cumI and bool_E_sum_I are True if there cumI and E+I are larger than 0, respectively
    ##                        <cumI>      <E+I>
    ##  not_yet_introduced      F           F
    ##  existing                T           T
    ##  extinction              T           F

    clade_not_introduced = (~bool_cumI)
    clade_existing = bool_cumI * bool_E_sum_I
    clade_extinction = bool_cumI * (~bool_E_sum_I)

    if data_flag == -1 : # before introduction of the clade (no data)
        if clade_not_introduced :
            return 1
        elif clade_existing :
            return 1
        else:
            return -1

    elif data_flag == 0 and np.isnan(data_S): # between the first and the last data points (either no data or yes data)
        if clade_not_introduced:
            return - 3
        elif clade_existing:
                return 1
        else:
            return -1


    elif data_flag == 0 and ~np.isnan(data_S): # between the first and the last data points (either no data or yes data)
        if clade_not_introduced:
            return - 3

        elif clade_existing:
            prob = poisson.pmf(data_S, particle_S)
            if np.isnan (prob):
                prob = -2
            return prob
        else:
            return -1

    elif data_flag == 1 : # after the last data points (no data)
        if clade_not_introduced:
            return -3
        elif clade_existing:
            return 1
        else:
            return 1





def get_weight_by_S_cmp (bool_cumI, bool_E_sum_I, data_flag, data_S, lambda_, nu_, logZ):

    # bool_cumI and bool_E_sum_I are True if there cumI and E+I are larger than 0, respectively
    ##                        <cumI>      <E+I>
    ##  not_yet_introduced      F           F
    ##  existing                T           T
    ##  extinction              T           F

    clade_not_introduced = (~bool_cumI)
    clade_existing = bool_cumI * bool_E_sum_I
    clade_extinction = bool_cumI * (~bool_E_sum_I)

    if data_flag == -1 : # before introduction of the clade (no data)
        if clade_not_introduced :
            return 1
        elif clade_existing :
            return 1
        else:
            return -1

    elif data_flag == 0 and np.isnan(data_S): # between the first and the last data points (either no data or yes data)
        if clade_not_introduced:
            return - 3
        elif clade_existing:
                return 1
        else:
            return -1


    elif data_flag == 0 and ~np.isnan(data_S): # between the first and the last data points (either no data or yes data)
        if clade_not_introduced:
            return - 3

        elif clade_existing:
            prob = cmp_pmf(data_S, lambda_, nu_, logZ)
            if np.isnan (prob):
                prob = -2
            return prob
        else:
            return -1

    elif data_flag == 1 : # after the last data points (no data)
        if clade_not_introduced:
            return -3
        elif clade_existing:
            return 1
        else:
            return 1



# -----------
# cmp related functions
def compute_logZ (log_lambda, nu, max_iter = 1e4, epsilon = 1e-8):
    # from https://github.com/jreduardo/cmpreg/blob/master/src/cmp_functions.cpp
    log_epsilon = np.log(epsilon)

    logz = 0
    logz_ = 0
    for j in range(1, int(max_iter)):
        logz_ += log_lambda - nu * np.log(j)
        logz = np.logaddexp(logz, logz_);

        if (logz_ - logz < log_epsilon) :
            break;

    return logz



def cmp_pmf (x, lambda_, nu_, logZ):

    #logZ = compute_logZ(np.log(lambda_), nu_)

    if x < 20:
        p = (lambda_ ** x) / ((factorial(x) ** nu_) * np.exp(logZ))

    if x < 100:
        log_p = x * np.log(lambda_) - nu_ * np.log(factorial(x)) - logZ
        p = np.exp(log_p)

    else:
        log_p = x * np.log(lambda_) - nu_ * (x * np.log(x) - x + 0.5 * np.log(2 * np.pi * x)) - logZ
        p = np.exp(log_p)

    return p

