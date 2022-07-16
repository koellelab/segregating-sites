import numpy as np
from numba import jit


def np_arange_prevent_imprecision(start, stop, step):
    # floating point overflow may cause the last element being greater than stop
    # -> To prevent this imprecision, we subtract (1e-4 * step) from the stop value

    return np.arange(start, stop - (step * 1e-10), step)




@jit(nopython = True)
def add_row_to_nparray(nparray, total_row_needed, increase_by = 1.5):

    result_nrow = nparray.shape[0]

    while True:
        result_nrow *= 1.5
        if int(result_nrow) > total_row_needed:
            break

    shape_to_add = (int(result_nrow) - nparray.shape[0] + 1, nparray.shape[1])

    #print ("shape_to_add:", shape_to_add)
    nparray = np.concatenate((nparray, np.zeros(shape_to_add, nparray.dtype)), axis = 0)

    return nparray




def np_arrange_with_endpoint(start, end, step):
  return np.arange(start, end + (step * 1e-10), step)