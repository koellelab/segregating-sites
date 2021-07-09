import numpy as np
from numba import jit


def np_arange_prevent_imprecision(start, stop, step):
    # floating point overflow may cause the last element being greater than stop
    # -> To prevent this imprecision, we subtract (1e-4 * step) from the stop value

    return np.arange(start, stop - (step * 1e-4), step)


def faster_random_sampling_without_replacement(element, num_samples, rng, probabilities = None, sample_size = 1):
    # https://medium.com/ibm-watson/incredibly-fast-random-sampling-in-python-baf154bd836a
    # a bit modified by yp(order of arguments, etc)

    '''
    :param element:       equivalent to <a> in np.random.choice
    :param num_samples:   equivalent to <size> in np.random.choice
    :param sample_size:   not in parameters for np.random.choice, but we can use this to get multiple sets of samples
    :param probabilities: equivalent to <p> in np.random.choice -> if uniform distirubution assumed, this will be np.ones(len(a))
                            todo - transmission heterogeneity??

    :return:
    '''

    if element.size == 1:
        element = np.arange(element)

    if probabilities == None:
        probabilities = np.ones_like(element) / element.size

    # this line is put inside of the function by yp for simplicity sake
    probabilities /= np.sum(probabilities)

    # replicate probabilities as many times as `num_samples`
    replicated_probabilities = np.tile(probabilities, (num_samples, 1))

    # get random shifting numbers & scale them correctly
    random_shifts = rng.random(replicated_probabilities.shape)
    random_shifts /= random_shifts.sum(axis=1)[:, np.newaxis]

    # shift by numbers & find largest (by finding the smallest of the negative)
    shifted_probabilities = random_shifts - replicated_probabilities

    return np.argpartition(shifted_probabilities, sample_size, axis=1)[:, :sample_size]


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


@jit(nopython=True)
def nparray_increase_row_col(nparray, total_row_needed, total_col_needed, increase_by = 1.5):

    arr_n_row = nparray.shape[0]
    arr_n_col = nparray.shape[1]

    while arr_n_row <= total_row_needed:
        arr_n_row *= increase_by

    while arr_n_col <= total_col_needed:
        arr_n_col *= increase_by

    nparray_new = np.zeros((int(arr_n_row)+1, int(arr_n_col)+1), nparray.dtype)
    nparray_new[:nparray.shape[0], :nparray.shape[1]] = nparray[:nparray.shape[0], :nparray.shape[1]]

    return nparray_new


@jit(nopython = True)
def tile_numba (array_1d, times):
    # https://stackoverflow.com/questions/61686293/numba-compatible-implementation-of-np-tile

    return np.repeat(array_1d, times).reshape(-1, times).T.flatten()

if __name__ == "__main__":

    ''' testing <faster_random_sampling_without_replacement>        '''
    # constants
    num_elements = 2000           # sampling from a population with 20 elements
    num_samples = 1          # 1000 sets of samples
    sample_size = 100             # each set of samples have 5
    elements = np.arange(num_elements)

    # probabilities should sum to 1
    count = np.arange(1, num_elements+1)
    probabilities = count / np.sum(count)

    sample = faster_random_sampling_without_replacement(num_samples, sample_size, probabilities)


    from matplotlib import pyplot as plt
    fig = plt.figure()
    plt.hist(sample.flatten()); #plt.xticks(np.arange(21))
    plt.show()



    pass
