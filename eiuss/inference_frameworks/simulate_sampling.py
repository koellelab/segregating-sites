import numpy as np
import pandas as pd

if __name__ == "__main__":
    import sys
    sys.path.append('../../eiuss')

import utils
import statistics_S as statistics

def _find_window(t_list, window):
    l = []
    for t in t_list:
        l.append(np.min(window[t <= window ], ))

    return np.array(l)

def _initiate_sample_data(data, win_intervals, n_clade=1):
    if not data is None:
        import copy
        sample_data = copy.deepcopy(data)
    else:
        sample_data = {}

    win_intervals_size = len(win_intervals)

    sample_data['win_intervals']      = win_intervals
    sample_data['n_samples_in_win']   = np.full([win_intervals_size, n_clade], np.nan)
    sample_data['n_recovered_in_win'] = np.full([win_intervals_size, n_clade], np.nan)
    sample_data['S_in_win']           = np.full([win_intervals_size, n_clade], np.nan)
    sample_data['sampled_individuals'] = np.empty([win_intervals_size, n_clade], dtype = 'object')

    return sample_data

def sampling_unif(sim, win_size, n_sample_per_window, rng,  sampling_duration = None):
    if sampling_duration:
        print(sampling_duration)
        sim_recovered_filter_by_time = (sampling_duration[0] < sim['recovered'][:, 0]) * (sim['recovered'][:, 0] <= sampling_duration[1])
        sim_recovered = sim['recovered'][sim_recovered_filter_by_time, :]
    else:
        sim_recovered = sim['recovered']
    sim_recovered[:, 0] = (sim_recovered[:, 0] * 1e6).astype('int') / 1e6

    state_var     = sim['statevar']
    win_intervals = utils.np_arange_prevent_imprecision(max(state_var[:, 0]), min(state_var[:, 0]) , win_size * (-1))[::-1]
    output_data   = _initiate_sample_data(None, win_intervals)

    ## find window for each recovered
    df = pd.DataFrame({'t': sim_recovered[:, 0],
                       'genotype': sim_recovered[:, 1],
                       'count': sim_recovered[:, 2],
                       'window': _find_window(sim_recovered[:, 0], win_intervals)})

    recovered_window = _find_window(sim_recovered[:, 0], win_intervals)


    clade_idx = 0
    sampled_df_cum = pd.DataFrame({'t': [], 'genotype': []})
    sampled_idx_cum = []
    for this_win_idx in range(len(win_intervals)):
        this_win_end = win_intervals[this_win_idx]
        recovered_idx_win = np.where(recovered_window == this_win_end)[0]
        recovered_idx_long = np.repeat(recovered_idx_win, sim_recovered[recovered_idx_win, 2].astype('int'))

        if len(recovered_idx_long) > n_sample_per_window:
            sampled_idx = rng.choice(recovered_idx_long, size=n_sample_per_window, replace=False)
            sampled_df = pd.DataFrame({'t': sim_recovered[sampled_idx, 0],
                                       'genotype': sim_recovered[sampled_idx, 1],
                                       'recovered_window': recovered_window[sampled_idx]})


            clade_genotype_info = sim['genotype_info'][:, 2 * clade_idx: 2 * (clade_idx + 1)]
            if len(sampled_df) > 1:
                output_data['S_in_win'][this_win_idx, clade_idx] = statistics.n_segregating_from_samples(
                    sampled_df['genotype'].astype('int'), clade_genotype_info)
                output_data['n_samples_in_win'][this_win_idx, clade_idx]   = len(sampled_df)
                output_data['n_recovered_in_win'][this_win_idx, clade_idx] = recovered_idx_long.size

                sampled_df_cum  = pd.concat([sampled_df_cum, sampled_df])
                sampled_idx_cum += sampled_idx.tolist()

    ## 4. save data
    output_data['sampled_idx']          = sampled_idx_cum
    output_data['sampled_individuals']  = sampled_df_cum
    output_data['win_intervals']        = win_intervals


    return output_data

def sampling_full(sim, win_size, rng):
    sim_recovered = sim['recovered']
    sim_recovered[:, 0] = (sim_recovered[:, 0] * 1e6).astype('int') / 1e6

    state_var = sim['statevar']
    win_intervals = utils.np_arange_prevent_imprecision(max(state_var[:, 0]), min(state_var[:, 0]) , win_size * (-1))[::-1]
    output_data   = _initiate_sample_data(None, win_intervals)

    ## find window for each recovered
    recovered_window = _find_window(sim_recovered[:, 0], win_intervals)

    clade_idx = 0
    for this_win_idx in range(len(win_intervals)):
        this_win_end = win_intervals[this_win_idx]
        recovered_idx_win = np.where(recovered_window == this_win_end)[0]
        recovered_idx_long = np.repeat(recovered_idx_win, sim['recovered'][recovered_idx_win, 2].astype('int'))

        if len(recovered_idx_long) > 0:
            sampled_idx = rng.choice(recovered_idx_long, size=len(recovered_idx_long), replace=False)
            sampled_df = pd.DataFrame({'t': sim['recovered'][sampled_idx, 0],
                                       'genotype': sim['recovered'][sampled_idx, 1],
                                      'recovered_window': recovered_window[sampled_idx]})

            clade_idx = 0
            clade_genotype_info = sim['genotype_info'][:, 2 * clade_idx: 2 * (clade_idx + 1)]
            if len(sampled_df) >= 1:
                output_data['S_in_win'][this_win_idx, clade_idx] = statistics.n_segregating_from_samples(
                    sampled_df['genotype'].astype('int'), clade_genotype_info)

                output_data['n_samples_in_win'][this_win_idx, clade_idx] = len(sampled_df)
                output_data['n_recovered_in_win'][this_win_idx, clade_idx] = recovered_idx_long.size


    return output_data

def sampling_prop(sim, win_size, n_total_sample, rng,  sampling_duration = None):
    if sampling_duration:
        print(sampling_duration)
        sim_recovered_filter_by_time = (sampling_duration[0] < sim['recovered'][:, 0]) * (sim['recovered'][:, 0] <= sampling_duration[1])
        sim_recovered = sim['recovered'][sim_recovered_filter_by_time, :]
    else:
        sim_recovered = sim['recovered']
    sim_recovered[:, 0] = (sim_recovered[:, 0] * 1e6).astype('int') / 1e6

    state_var     = sim['statevar']
    win_intervals = utils.np_arange_prevent_imprecision(max(state_var[:, 0]) + win_size*0.5, min(state_var[:, 0]) , win_size * (-1))[::-1]
    output_data   = _initiate_sample_data(None, win_intervals)

    ## find window for each recovered
    recovered_window = _find_window(sim_recovered[:, 0], win_intervals)
    df = pd.DataFrame({'t': sim_recovered[:, 0],
                       'genotype': sim_recovered[:, 1],
                       'count': sim_recovered[:, 2],
                       'window': _find_window(sim_recovered[:, 0], win_intervals)})

    ## get n_samples for each window
    window_for_ind = np.repeat(recovered_window, sim_recovered[:, 2].astype('int'))
    selected_win = rng.choice(window_for_ind, size = n_total_sample, replace=False)
    window_end, n_sample_for_window = np.unique(selected_win, return_counts=True)

    for i in range(len(n_sample_for_window)):
        print(window_end[i], n_sample_for_window[i])


    clade_idx = 0
    sampled_df_cum = pd.DataFrame({'t':[], 'genotype':[]})
    sampled_idx_cum = []
    for this_win_idx in range(len(win_intervals)):
        this_win_end = win_intervals[this_win_idx]
        recovered_idx_win = np.where(recovered_window == this_win_end)[0]
        recovered_idx_long = np.repeat(recovered_idx_win, sim_recovered[recovered_idx_win, 2].astype('int'))

        ## get n_sample
        n_count_idx = np.where(window_end == this_win_end)[0]
        if n_count_idx.size == 0:
            n_sample_per_window = 0
        else:
            n_sample_per_window = n_sample_for_window[n_count_idx][0]



        if len(recovered_idx_long) >= n_sample_per_window:
            sampled_idx = rng.choice(recovered_idx_long, size=n_sample_per_window, replace=False)
            sampled_df = pd.DataFrame({'t': sim_recovered[sampled_idx, 0],
                                       'genotype': sim_recovered[sampled_idx, 1],
                                       'recovered_window': recovered_window[sampled_idx]})

            clade_genotype_info = sim['genotype_info'][:, 2 * clade_idx: 2 * (clade_idx + 1)]
            if len(sampled_df) >= 1:
                #print(f"{this_win_end:5.2f}, {n_sample_per_window:7d}, {len(recovered_idx_long):7d} v ")
                output_data['S_in_win'][this_win_idx, clade_idx] = statistics.n_segregating_from_samples(
                    sampled_df['genotype'].astype('int'), clade_genotype_info)
                output_data['n_samples_in_win'][this_win_idx, clade_idx] = len(sampled_df)
                output_data['n_recovered_in_win'][this_win_idx, clade_idx] = recovered_idx_long.size

                sampled_df_cum = pd.concat([sampled_df_cum, sampled_df])
                sampled_idx_cum += sampled_idx.tolist()

        else:
            print(f"{this_win_end:5.2f}, {n_sample_per_window:7d}, {len(recovered_idx_long):7d}  ")
        ## 4. save data
        output_data['sampled_idx'] = sampled_idx_cum
        output_data['sampled_individuals'] = sampled_df_cum
        output_data['win_intervals'] = win_intervals




    return output_data


def sampling_prop_large_set(sim, win_size, n_total_sample, rng,  sampling_duration = None):
    if sampling_duration:
        print(sampling_duration)
        sim_recovered_filter_by_time = (sampling_duration[0] < sim['recovered'][:, 0]) * (sim['recovered'][:, 0] <= sampling_duration[1])
        sim_recovered = sim['recovered'][sim_recovered_filter_by_time, :]
    else:
        sim_recovered = sim['recovered']
    sim_recovered[:, 0] = (sim_recovered[:, 0] * 1e6).astype('int') / 1e6


    recovered_idx_long = np.repeat(np.arange(sim_recovered[:, 2].size), sim_recovered[:, 2].astype('int'))
    n_total_recovered = len(recovered_idx_long)
    print(f"n_total_recovered = {n_total_recovered}")

    # now sampling, n_total_sample in each clade
    for clade_idx in range(sim['n_clade']):
        clade_genotype_info = sim['genotype_info'][:, 2 * clade_idx: 2 * (clade_idx + 1)]

        ## 1. sample recovered individuals
        sampled_idx = rng.choice(recovered_idx_long, size=n_total_sample, replace=False)
        sampled_df = pd.DataFrame({'t'       : sim_recovered[sampled_idx, 0],
                                   'genotype': sim_recovered[sampled_idx, 1]}).reset_index(drop=True)

        ## 2. redefine the window interval using the last sampling time
        last_sample_time = sampled_df['t'].max()
        win_intervals = utils.np_arange_prevent_imprecision(        ## window interval is redefined with last sampling time
            np.ceil(last_sample_time),
            sim['statevar'][:, 0].min(),
            win_size * (-1))[::-1]

        ## 3. find which window each sample came from
        print("last sample time")
        sampled_df['window'] = win_intervals[np.digitize(sampled_df['t'],     bins=win_intervals, right = True)]
        ## also, count the number of recovered in each window
        recovered_time          = sim['recovered'][recovered_idx_long, 0]
        recovered_time_filtered = recovered_time[recovered_time <= np.ceil(last_sample_time)]
        recovered_window_window = win_intervals[np.digitize(recovered_time_filtered, bins=win_intervals, right=True)]
        #print(f"see what is not sampled: {recovered_time[recovered_time > last_sample_time]}")
        sim.close()

        ## 4. save data
        output_data = _initiate_sample_data(None, win_intervals)
        output_data['sampled_idx']         = sampled_idx
        output_data['sampled_individuals'] = sampled_df
        output_data['win_intervals']       = win_intervals

        ## 5. for each window, calculate the statistics
        for this_win_idx in range(len(win_intervals)):
            print("this_win_idx:", this_win_idx, "/", len(win_intervals))
            this_win_end = win_intervals[this_win_idx]

            this_win_sample     = sampled_df.loc[sampled_df.loc[:, 'window'] == this_win_end]
            this_win_recovered  = (recovered_window_window == this_win_end).sum()
            if len(this_win_sample) > 0:
                output_data['S_in_win'][this_win_idx, clade_idx] = statistics.n_segregating_from_samples(
                    this_win_sample['genotype'].to_numpy().astype('int'), clade_genotype_info)
            output_data['n_samples_in_win'][this_win_idx, clade_idx] = len(this_win_sample)
            output_data['n_recovered_in_win'][this_win_idx, clade_idx] = this_win_recovered

    return output_data





def sampling_reset_window(sim, win_size, input_samp, rng,  sampling_duration = None):
    samp = np.load(f'{input_samp}_sampling.npz', allow_pickle=True)
    #samp = np.load("/Users/yeongseon/Dropbox/PhD_Emory_university/1_projects/Project_seg_site/simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_sampling.npz", allow_pickle=True)

    if sampling_duration:
        print(sampling_duration)
        sim_recovered_filter_by_time = (sampling_duration[0] < sim['recovered'][:, 0]) * (sim['recovered'][:, 0] <= sampling_duration[1])
        sample_filter_by_time = (sampling_duration[0] < samp['sampled_individuals'][:, 0]) * (samp['sampled_individuals'][:, 0] <= sampling_duration[1])

        sim_recovered = sim['recovered'][sim_recovered_filter_by_time, :]
        already_sampled = samp['sampled_individuals'][sample_filter_by_time, :]
    else:
        sim_recovered = sim['recovered']
        already_sampled = samp['sampled_individuals']

    sim_recovered[:, 0] = (sim_recovered[:, 0] * 1e6).astype('int') / 1e6



    state_var     = sim['statevar']
    win_intervals = utils.np_arange_prevent_imprecision(max(state_var[:, 0]) + win_size*0.5, min(state_var[:, 0]) , win_size * (-1))[::-1]
    output_data   = _initiate_sample_data(None, win_intervals)

    ## get new recovered_window for already sampled individuals
    sampled_window   = _find_window(already_sampled[:, 0], win_intervals)
    recovered_window = _find_window(  sim_recovered[:, 0], win_intervals)

    clade_idx = 0
    sampled_df_cum = pd.DataFrame({'t':[], 'genotype':[]})
    sampled_idx_cum = []
    for this_win_idx in range(len(win_intervals)):
        this_win_end = win_intervals[this_win_idx]


        recovered_idx_win = np.where(recovered_window == this_win_end)[0]
        recovered_idx_long = np.repeat(recovered_idx_win, sim_recovered[recovered_idx_win, 2].astype('int'))


        sampled_idx = np.where(sampled_window == this_win_end)[0]
        sampled_df = pd.DataFrame(already_sampled[sampled_idx,:2], columns = ['t', 'genotype'])
        sampled_df['recovered_window'] = sampled_window[sampled_idx]
        print(sampled_df.to_string())

        ## get segregating sites
        clade_genotype_info = sim['genotype_info'][:, 2 * clade_idx: 2 * (clade_idx + 1)]
        if len(sampled_df) >= 1:
            output_data['S_in_win'][this_win_idx, clade_idx] = statistics.n_segregating_from_samples(
                sampled_df['genotype'].astype('int'), clade_genotype_info)
            output_data['n_samples_in_win'][this_win_idx, clade_idx] = len(sampled_df)
            output_data['n_recovered_in_win'][this_win_idx, clade_idx] = recovered_idx_long.size #sim_recovered[recovered_window == this_win_end, 2].sum()

            sampled_df_cum = pd.concat([sampled_df_cum, sampled_df])
            sampled_idx_cum += sampled_idx.tolist()


    ## 4. save data
    output_data['sampled_idx'] = sampled_idx_cum
    output_data['sampled_individuals'] = sampled_df_cum
    output_data['win_intervals'] = win_intervals


    return output_data









def run (params, none1, none2):
    sim = np.load(f'{params["input"]}_simtaj.npz', allow_pickle=True)


    if params['sample_type'] == "unif":
        #data = sampling_unif(sim, params['window_dt'], params['sample_size'] , np.random.default_rng(params['seed']), params['sample_between'])
        data = sampling_unif(sim, params['window_dt'], params['sample_size'], params['rng'], params['sample_between'])

    elif params['sample_type'] == "full":
        #data = sampling_full(sim, params['window_dt'], np.random.default_rng(params['seed']))
        data = sampling_full(sim, params['window_dt'],  params['rng'])

    elif params['sample_type'] == "prop":
        #data = sampling_prop(sim, params['window_dt'], params['sample_size'], np.random.default_rng(params['seed']), params['sample_between'])
        data = sampling_prop(sim, params['window_dt'], params['sample_size'], params['rng'],params['sample_between'])

    elif params['sample_type'] == "reset_window":
        data = sampling_reset_window(sim, params['window_dt'], params['input_samp'], params['rng'], params['sample_between'])










    np.savez(params['out_name'] + "_sampling.npz", **{key: value for key, value in data.items()})
    dir_table = params['out_name'] + "_segsites.tsv"
    f = open(dir_table, 'w')
    print('\t'.join(["window", "n", "s", "window_end", "n_recovered_in_win", "clade"]), file=f)

    n_clade = 1
    n_sample_sum = [0] * n_clade
    for i in range(len(data['win_intervals'])):
        for clade_idx in range(n_clade):
            data_values = [i, data['n_samples_in_win'][i][clade_idx], data['S_in_win'][i][clade_idx],  data['win_intervals'][i], data['n_recovered_in_win'][i][clade_idx], clade_idx]
            print('\t'.join([str(x) for x in data_values]))


            if ~np.isnan(data['n_samples_in_win'][i][clade_idx]):
                print('\t'.join([str(x) for x in data_values]), file=f, flush=True)
                n_sample_sum[clade_idx] += data['n_samples_in_win'][i][clade_idx]

    print ("n_sample_sum: ", n_sample_sum)
    data['sampled_individuals'].to_csv(params['out_name'] + "_sampled_individuals.tsv", sep = "\t")


    return



if __name__ == "__main__":


    params = {}
    params['sample_type'] = 'reset_window'
    params['window_dt'] = 4
    params['sample_size'] = 500
    params['rng']   =  np.random.default_rng(12345)
    params['sample_between'] = None


    params['input']      = "/Users/yeongseon/Dropbox/PhD_Emory_university/1_projects/Project_seg_site/simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234"
    params['input_samp'] = "/Users/yeongseon/Dropbox/PhD_Emory_university/1_projects/Project_seg_site/simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201"
    params['out_name'] =params['input'] + "_test"
    run (params, None, None)

