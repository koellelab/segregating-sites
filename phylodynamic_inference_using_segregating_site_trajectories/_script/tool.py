# ------------------
## tool.py 
## Author : Yeongseon Park
##          Mike Martin 
# ------------------
import numpy as np
import pandas as pd
from numba import jit
import scipy.io as spio
from datetime import datetime
from matplotlib import pyplot as plt


def np_arrange_with_endpoint(start, end, step):
  return np.arange(start, end + (step * 0.001), step)



def open_benchmark_log(filename_benchmark, write_mode='a'):
    benchmark_file = open(filename_benchmark, write_mode)

    if write_mode == "w":
        print("* Please check whether random number generators are controlled to get the accurate benchmark result", file=benchmark_file)
        print("\n\n\n" + "=" * 146, file=benchmark_file)

        benchmark_table = ["repeats", "max", "min", "mean", "std", "median", "version_nickname"]
        benchmark_table = ["{0:>24s}".format("time")] + ["{0:>12s}".format(x) for x in benchmark_table]
        benchmark_table = "\t".join(benchmark_table)

        print(benchmark_table, file=benchmark_file)
        print("-" * 146, file=benchmark_file)

    return benchmark_file


def write_benchmark_log(file_benchmark, time_record, version_nickname="-", save_data=False):
    benchmark_result = [time_record.size, time_record.max(), time_record.min(), time_record.mean(), time_record.std(), np.median(time_record)]
    benchmark_result = ["{0: >#12.6f}".format(x) for x in benchmark_result]
    benchmark_result = datetime.now().ctime() + "\t" + "\t".join(benchmark_result) + "\t" + version_nickname

    print(benchmark_result, file=file_benchmark, flush=True)

    if save_data:
        # filename_data = open(file_benchmark.name + version_nickname, "w")
        filename_data = "/Users/yeongseon/bt.txt"

        try:
            file_data = open(filename_data, "x")
        except FileExistsError:  # if there is already a folder start with timestamp
            print("[Benchmark ERROR]: version_nickname already exists. Timestamp is used instead")
            filename_data = filename_data + "." + datetime.now().strftime('%y%m%d-%H%M$S')
            file_data = open(filename_data, "x")
        finally:
            print(time_record, sep="\n", file=file_data, flush=True)

    return




@jit(nopython = True)
def tile_numba (array_1d, times):
    # https://stackoverflow.com/questions/61686293/numba-compatible-implementation-of-np-tile

    return np.repeat(array_1d, times).reshape(-1, times).T.flatten()


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


def get_figure_SEIR_day(particle, params, n_clade):
    (fig, axes) = plt.subplots(2, 2, figsize=[10, 8])
    state_var = particle.statevar[particle.statevar[:, 0] > 0, :]


    axes[0, 0].plot(state_var[:, 0], state_var[:, 1])
    for i in range(n_clade):
        print (i)
        axes[0, 1].plot(state_var[:, 0], state_var[:, i * 3 + 2], label="clade #" + str(i))
        axes[1, 0].plot(state_var[:, 0], state_var[:, i * 3 + 3], label="clade #" + str(i))
        axes[1, 1].plot(state_var[:, 0], state_var[:, i * 3 + 4], label="clade #" + str(i))

    axes[0, 0].set_ylabel('S')
    axes[0, 1].set_ylabel('E')
    axes[1, 0].set_ylabel('I')
    axes[1, 1].set_ylabel('cumI')

    axes[0, 0].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])
    axes[0, 1].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])
    axes[1, 0].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])
    axes[1, 1].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])

    axes[0, 0].set_ylim([0, params['N']])
    axes[1, 1].set_ylim([0, params['N']])

    axes[1, 1].legend();
    fig.tight_layout();

    return fig


def import_data(data_path):
    # data -> copied from Mike's script

    dat_dict = {}
    dat = pd.read_csv(data_path, sep='\t')

    # these items are the same for each clade
    # assumes windows are evenly spaced
    # assumes all clades have the same windows
    len_win = len(dat.iloc[:, :]['window_end']); print ("len_win:", len_win)
    dat_dict['window_dt'] = np.min(np.array(dat.iloc[1:len_win, :]['window_end']) - np.array(dat.iloc[0:len_win-1, :]['window_end']))
    print ("window_dt:", dat_dict['window_dt'])
    dat_dict['win_intervals'] = np.array(sorted(set(dat['window_end'])))
    clades = set(dat['clade'])
    dat_dict['n_clades'] = len(clades)

    # defines empty arrays
    dat_dict['nsequences'] = np.empty(shape=(len(clades), 1))
    dat_dict['taj_d_in_win'] = np.empty(shape=(len(dat_dict['win_intervals']), len(clades)))
    dat_dict['n_samples_in_win'] = np.empty(shape=(len(dat_dict['win_intervals']), len(clades)))
    dat_dict['pi_in_win'] = np.empty(shape=(len(dat_dict['win_intervals']), len(clades)))
    dat_dict['S_in_win'] = np.empty(shape=(len(dat_dict['win_intervals']), len(clades)))
    dat_dict['var_taj_d_in_win'] = np.empty(shape=(len(dat_dict['win_intervals']), len(clades)))        ## added by yp

    # fiills arrays with clade data
    for clade_idx, clade in enumerate(clades):
        clade_dat = dat[dat['clade'] == clade]
        for t_idx, t in enumerate(dat_dict['win_intervals']):
            dat_dict['taj_d_in_win'][t_idx, clade_idx] = clade_dat['taj_d'][t_idx]
            dat_dict['n_samples_in_win'][t_idx, clade_idx] = clade_dat['n'][t_idx]
            dat_dict['pi_in_win'][t_idx, clade_idx] = clade_dat['pi'][t_idx]
            dat_dict['S_in_win'][t_idx, clade_idx] = clade_dat['s'][t_idx]

            if dat_dict['n_samples_in_win'][t_idx, clade_idx] > 2 :         ## added by yp
                dat_dict['var_taj_d_in_win'][t_idx, clade_idx] = 1          ## added by yp


    return (dat_dict)







def import_data(data_path):
    # data -> copied from Mike's script

    dat_dict = {}
    dat = pd.read_csv(data_path, sep='\t')

    # these items are the same for each clade
    # assumes windows are evenly spaced
    # assumes all clades have the same windows
    len_win = len(dat.iloc[:, :]['window_end']); print ("len_win:", len_win)
    dat_dict['window_dt'] = np.min(np.array(dat.iloc[1:len_win, :]['window_end']) - np.array(dat.iloc[0:len_win-1, :]['window_end']))
    print ("window_dt:", dat_dict['window_dt'])
    dat_dict['win_intervals'] = np.array(sorted(set(dat['window_end'])))
    clades = set(dat['clade'])
    dat_dict['n_clades'] = len(clades)

    # defines empty arrays
    dat_dict['nsequences'] = np.empty(shape=(len(clades), 1))
    dat_dict['taj_d_in_win'] = np.empty(shape=(len(dat_dict['win_intervals']), len(clades)))
    dat_dict['n_samples_in_win'] = np.empty(shape=(len(dat_dict['win_intervals']), len(clades)))
    dat_dict['pi_in_win'] = np.empty(shape=(len(dat_dict['win_intervals']), len(clades)))
    dat_dict['S_in_win'] = np.empty(shape=(len(dat_dict['win_intervals']), len(clades)))
    dat_dict['var_taj_d_in_win'] = np.empty(shape=(len(dat_dict['win_intervals']), len(clades)))        ## added by yp

    # fiills arrays with clade data
    for clade_idx, clade in enumerate(clades):
        clade_dat = dat[dat['clade'] == clade]
        for t_idx, t in enumerate(dat_dict['win_intervals']):
            dat_dict['taj_d_in_win'][t_idx, clade_idx] = clade_dat['taj_d'][t_idx]
            dat_dict['n_samples_in_win'][t_idx, clade_idx] = clade_dat['n'][t_idx]
            dat_dict['pi_in_win'][t_idx, clade_idx] = clade_dat['pi'][t_idx]
            dat_dict['S_in_win'][t_idx, clade_idx] = clade_dat['s'][t_idx]

            if dat_dict['n_samples_in_win'][t_idx, clade_idx] > 2 :         ## added by yp
                dat_dict['var_taj_d_in_win'][t_idx, clade_idx] = 1          ## added by yp


    return (dat_dict)




# from https://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict


def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict







def import_data_n_segregating(data_path):
    # data -> copied from Mike's script

    dat_dict = {}
    dat = pd.read_csv(data_path, sep='\t')

    # these items are the same for each clade
    # assumes windows are evenly spaced
    # assumes all clades have the same windows
    len_win = len(dat.iloc[:, :]['window_end']); print ("len_win:", len_win)
    dat_dict['window_dt'] = np.min(np.array(dat.iloc[1:len_win, :]['window_end']) - np.array(dat.iloc[0:len_win-1, :]['window_end']))
    print ("window_dt:", dat_dict['window_dt'])
    dat_dict['win_intervals'] = np.array(sorted(set(dat['window_end'])))
    clades = set(dat['clade'])
    dat_dict['n_clades'] = len(clades)

    # defines empty arrays
    dat_dict['nsequences'] = np.empty(shape=(len(clades), 1))
    dat_dict['n_samples_in_win'] = np.empty(shape=(len(dat_dict['win_intervals']), len(clades)))
    dat_dict['S_in_win'] = np.empty(shape=(len(dat_dict['win_intervals']), len(clades)))

    # fiills arrays with clade data
    for clade_idx, clade in enumerate(clades):
        clade_dat = dat[dat['clade'] == clade]
        for t_idx, t in enumerate(dat_dict['win_intervals']):

            dat_dict['n_samples_in_win'][t_idx, clade_idx] = clade_dat['n'][t_idx]
            dat_dict['S_in_win'][t_idx, clade_idx] = clade_dat['s'][t_idx]




    return (dat_dict)






if __name__ == "__main__":
    a = np.ones([3, 6])

    b = nparray_increase_row_col(a, 5, 18)
    print (b)