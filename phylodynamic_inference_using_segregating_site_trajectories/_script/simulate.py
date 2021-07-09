import time
import numpy as np
import pandas as pd
from tools import misc_
from tools import figure

import statistics_S as statistics


def SE1E2IR_hetero_exact_Salje (params_sim, dataset_dir, file_header):

    from SE1E2IR_hetero import next_dt, particle

    sim = particle.Particle(n_clade=1, n_class = params_sim['n_class'], init_size= 1201)
    sim_window = misc_.np_arange_prevent_imprecision(params_sim['timeend'], params_sim['timeintro'], -1 * params_sim['dt'])[::-1]



    while True:
        time_start = time.time()

        sim.initiate(params_sim, 1201)
        sim.recovered = pd.DataFrame()
        sim.transhet = np.zeros(5, dtype = "float")

        # ------------------------------------------------
        # in Salje, initial cases start with E = 58.65
        sim.statevar[0, 1] -= 57            # S
        sim.statevar[0, 2] += 57            # E
        sim.statevar[0, 7] += 57            # cumI

        sim.genotype_count[0, 0] += 57      # E
        # ------------------------------------------------

        this_time_start = params_sim['timeintro']
        for idx in range(len(sim_window)):      # = 1200
            this_time_end = sim_window[idx]

            #next_dt.update (sim, params_sim, this_time_start, this_time_end)
            next_dt.update_track_transhet(sim, params_sim, this_time_start, this_time_end)

            for clade_idx in range(sim.n_clade):
                n_genotype_max = np.max(sim.n_genotype)
                idx_recovered = 2 + 3 * clade_idx

                #print(add_df)

            this_time_start = this_time_end


        # transmission heterogeneity
        high_group = sim.transhet[2] / (sim.transhet[1] + sim.transhet[2])              # E1_to_E2h / (E1_to_E2h + E1_to_E2l)
        infeced_by_high = sim.transhet[4] / (sim.transhet[3] + sim.transhet[4])         # S_to_E_Ih / (S_to_E_Ih + S_to_E_Il)



        printStr = "n_genotype: " + str(sim.n_genotype) + ", n_site: " + str(sim.n_site) + "\n"
        printStr +=  "%.3f (%d) of high are responsible for \n %.3f (%d) of secondary infection\n" %(high_group, sim.transhet[2], infeced_by_high, sim.transhet[4])
        print(printStr)




        print("runtime: ", time.time() - time_start, "(sec)")
        sim.trim_arrays()

        fig1 = figure.get_figure_SE1E2IR_day_oneclade(sim, params_sim, printStr)
        fig1.savefig(dataset_dir + file_header + ".png")

        fig2 = figure.get_figure_SE1E2IR_day_oneclade_type2(sim, params_sim, printStr)
        fig2.savefig(dataset_dir + file_header + ".png")

        print(sim.recovered)
        np.savez(dataset_dir + file_header + "_simtaj.npz", **{key: np.array(value) for key, value in sim.__dict__.items()})


        break





    return sim




def SE1E2IR_hetero (params_sim, dataset_dir, file_header):

    from SE1E2IR_hetero import next_dt, particle

    sim = particle.Particle(n_clade=1, n_class = params_sim['n_class'], init_size= 1201)
    sim_window = misc_.np_arange_prevent_imprecision(params_sim['timeend'], params_sim['timeintro'], -1 * params_sim['dt'])[::-1]

    while True:
        time_start = time.time()

        sim.initiate(params_sim, 1201)
        sim.recovered = pd.DataFrame()
        sim.transhet = np.zeros(5, dtype = "float")

        this_time_start = params_sim['timeintro']
        for idx in range(len(sim_window)):      # = 1200
            this_time_end = sim_window[idx]

            #next_dt.update (sim, params_sim, this_time_start, this_time_end)
            next_dt.update_track_transhet(sim, params_sim, this_time_start, this_time_end)

            for clade_idx in range(sim.n_clade):
                n_genotype_max = np.max(sim.n_genotype)
                idx_recovered = 2 + 3 * clade_idx

                idx = sim.genotype_count[:n_genotype_max, idx_recovered] > 0
                sim.recovered = sim.recovered.append(pd.DataFrame({'t': this_time_end,
                                                                   'genotype': np.arange(n_genotype_max)[idx],
                                                                   'count': sim.genotype_count[:n_genotype_max, idx_recovered][idx],
                                                                   'clade': clade_idx}))

                #print(add_df)

            this_time_start = this_time_end


        # transmission heterogeneity
        high_group = sim.transhet[2] / (sim.transhet[1] + sim.transhet[2])              # E1_to_E2h / (E1_to_E2h + E1_to_E2l)
        infeced_by_high = sim.transhet[4] / (sim.transhet[3] + sim.transhet[4])         # S_to_E_Ih / (S_to_E_Ih + S_to_E_Il)



        printStr = "n_genotype: " + str(sim.n_genotype) + ", n_site: " + str(sim.n_site) + "\n"
        printStr +=  "%.3f of high are responsible for \n %.3f of secondary infection\n" %(high_group, infeced_by_high)
        print(printStr)



        if (100 < sim.n_genotype).all():  # and (sim.n_genotype < 1000).all():
            print("runtime: ", time.time() - time_start, "(sec)")
            sim.trim_arrays()

            fig1 = figure.get_figure_SE1E2IR_day_oneclade(sim, params_sim, printStr)
            fig1.savefig(dataset_dir + file_header + ".png")

            fig2 = figure.get_figure_SE1E2IR_day_oneclade_type2(sim, params_sim, printStr)
            fig2.savefig(dataset_dir + file_header + ".png")

            print(sim.recovered)
            sim.recovered = sim.recovered.to_numpy()

            np.savez(dataset_dir + file_header + "_simtaj.npz",**{key: np.array(value) for key, value in sim.__dict__.items()})
            np.savez(dataset_dir + file_header + "_params.npz",**{key: np.array(value) for key, value in params_sim.items()})

            break


    return



def SEIR_simple (params_sim, dataset_dir, file_header):

    from SEIR_simple import next_dt, particle

    sim = particle.Particle(n_clade=1, n_class = params_sim['n_class'], init_size= 1201)
    sim_window = misc_.np_arange_prevent_imprecision(params_sim['timeend'], params_sim['timeintro'], -1 * params_sim['dt'])[::-1]

    while True:
        time_start = time.time()

        sim.initiate(params_sim, 1201)
        sim.recovered = pd.DataFrame()

        this_time_start = params_sim['timeintro']
        for idx in range(len(sim_window)):      # = 1200
            this_time_end = sim_window[idx]

            #next_dt.update (sim, params_sim, this_time_start, this_time_end)
            next_dt.update(sim, params_sim, this_time_start, this_time_end)

            for clade_idx in range(sim.n_clade):
                n_genotype_max = np.max(sim.n_genotype)
                idx_recovered = 2 + 3 * clade_idx

                idx = sim.genotype_count[:n_genotype_max, idx_recovered] > 0
                sim.recovered = sim.recovered.append(pd.DataFrame({'t': this_time_end,
                                                                   'genotype': np.arange(n_genotype_max)[idx].astype('int'),
                                                                   'count': sim.genotype_count[:n_genotype_max, idx_recovered][idx].astype('int'),
                                                                   'clade': clade_idx}))

                #print(add_df)

            this_time_start = this_time_end

        printStr = "n_genotype: " + str(sim.n_genotype) + ", n_site: " + str(sim.n_site) + "\n"
        print(printStr)

        if (100 < sim.n_genotype).all():  # and (sim.n_genotype < 1000).all():
            print("runtime: ", time.time() - time_start, "(sec)")
            sim.trim_arrays()

            fig1 = figure.get_figure_SEIR_day_oneclade(sim, params_sim, printStr)
            fig1.savefig(dataset_dir + file_header + ".png")

            print(sim.recovered)
            sim.recovered = sim.recovered.to_numpy()

            np.savez(dataset_dir + file_header + "_simtaj.npz", **{key: np.array(value) for key, value in sim.__dict__.items()})
            np.savez(dataset_dir + file_header + "_params.npz", **{key: np.array(value) for key, value in params_sim.items()})

            break


    return





import numpy as np
import pandas as pd

def sampling_prop (sim, data, n_total_sample):

    win_intervals = data['win_intervals']
    data = _initiate_sample_data(data, len(win_intervals), sim['n_clade'])



    recovered_df = pd.DataFrame({'t': sim['recovered'][:, 0],
                                 'genotype': sim['recovered'][:, 1],
                                 'count': sim['recovered'][:, 2],
                                 'clade': sim['recovered'][:, 3],
                                 'window': _find_window(sim['recovered'][:, 0], win_intervals)})


    recovered_idx_long = np.repeat(np.arange(len(recovered_df)), recovered_df['count'])


    #recovered_df = recovered_df.iloc[np.repeat(np.arange(len(recovered_df)), recovered_df['count'])].reset_index(drop = True)


    # now sampling, n_total_sample in each clade
    for clade_idx in range(sim['n_clade']):

        clade_genotype_info = sim['genotype_info'][:, 2 * clade_idx: 2 * (clade_idx + 1)]

        #sampled_idx = np.random.choice(len(recovered_df), size = n_total_sample)
        sampled_idx = np.random.choice(recovered_idx_long, size=n_total_sample)
        sampled_df = recovered_df.iloc[sampled_idx].reset_index(drop = True)

        for this_win_idx in range(len(win_intervals)):
            print ("this_win_idx:", this_win_idx, "/", len(win_intervals))
            this_win_end = win_intervals[this_win_idx]

            this_win_sample = sampled_df.loc[sampled_df.loc[:, 'window'] == this_win_end]

            if len(this_win_sample) > 1:
                data['S_in_win'][this_win_idx, clade_idx] = statistics.n_segregating_from_samples(this_win_sample['genotype'].astype('int'), clade_genotype_info)

            data['n_samples_in_win'][this_win_idx, clade_idx] = len(this_win_sample)

    return data


def sampling_prop_large_set(sim, data, n_total_sample):
    win_intervals = data['win_intervals']
    data = _initiate_sample_data(data, len(win_intervals), sim['n_clade'])

    # single clade approach
    recovered_idx_long = np.repeat(np.arange(sim['recovered'][:, 2].size), sim['recovered'][:, 2].astype('int'))
    # now sampling, n_total_sample in each clade
    for clade_idx in range(sim['n_clade']):

        clade_genotype_info = sim['genotype_info'][:, 2 * clade_idx: 2 * (clade_idx + 1)]

        sampled_idx = np.random.choice(recovered_idx_long, size=n_total_sample)
        sampled_df = pd.DataFrame({'t': sim['recovered'][sampled_idx, 0],
                                   'genotype': sim['recovered'][sampled_idx, 1],
                                   'count': sim['recovered'][sampled_idx, 2],
                                   'clade': sim['recovered'][sampled_idx, 3]})
        sim.close()
        sampled_df['window'] = _find_window(sampled_df['t'], win_intervals)


        for this_win_idx in range(len(win_intervals)):
            print("this_win_idx:", this_win_idx, "/", len(win_intervals))
            this_win_end = win_intervals[this_win_idx]

            this_win_sample = sampled_df.loc[sampled_df.loc[:, 'window'] == this_win_end]

            #if len(this_win_sample) > 0:
            #    data['S_in_win'][this_win_idx, clade_idx] = statistics.n_segregating_from_samples(
            #        this_win_sample['genotype'].to_numpy.astype('int'), clade_genotype_info)

            data['n_samples_in_win'][this_win_idx, clade_idx] = len(this_win_sample)

    return data


def sampling_data_large_set (sim, data, imported_data):
    win_intervals = data['win_intervals']
    data = _initiate_sample_data(data, len(win_intervals), sim['n_clade'])

    recovered_window = _find_window(sim['recovered'][:, 0], win_intervals)

    for this_win_idx in range(len(win_intervals)):
        print("this_win_idx:", this_win_idx, "/", len(win_intervals))
        this_win_end = win_intervals[this_win_idx]

        this_win_n_sample = imported_data.loc[imported_data['win_end'] == this_win_end]['n_sample'].values.astype('int')
        print (this_win_n_sample)

        if this_win_n_sample > 0:

            recovered_idx_win = np.where(recovered_window == this_win_end)[0]
            recovered_idx_long = np.repeat(recovered_idx_win, sim['recovered'][recovered_idx_win, 2].astype('int'))
            sampled_idx = np.random.choice(recovered_idx_long, size=this_win_n_sample)

            sampled_df = pd.DataFrame({'t': sim['recovered'][sampled_idx, 0],
                                         'genotype': sim['recovered'][sampled_idx, 1],
                                         'count': sim['recovered'][sampled_idx, 2],
                                         'clade': sim['recovered'][sampled_idx, 3]})

            clade_idx = 0
            clade_genotype_info = sim['genotype_info'][:, 2 * clade_idx: 2 * (clade_idx + 1)]
            if len(sampled_df) > 1:
                data['S_in_win'][this_win_idx, clade_idx] = statistics.n_segregating_from_samples(
                    sampled_df['genotype'].astype('int'), clade_genotype_info)

                data['n_samples_in_win'][this_win_idx, clade_idx] = this_win_n_sample

    return data


def sampling_prop_with_max(sim, data, n_total_sample, max_sample_size):
    win_intervals = data['win_intervals']
    data = _initiate_sample_data(data, len(win_intervals), sim['n_clade'])

    recovered_window = _find_window(sim['recovered'][:, 0], win_intervals)

    # 1. get the number of recovered individuals in each window
    recovered_idx = []
    recovered_count = []
    for this_win_idx in range(len(win_intervals)):
        this_win_end = win_intervals[this_win_idx]
        recovered_idx_win = np.where(recovered_window == this_win_end)[0]

        recovered_idx.append(recovered_idx_win.tolist())
        recovered_count.append(recovered_idx_win.size)

    # 2. get the number of samples in each window
    recovered_count_orig = np.array(recovered_count)
    recovered_count = np.array(recovered_count)
    recovered_freq = recovered_count / recovered_count.sum()
    print("recovered_count:", recovered_count)
    print("recovered_freq:", recovered_freq)

    n_remaining_sample = n_total_sample
    n_sample_win = np.full(len(win_intervals), 0, dtype='int')

    ## the number of maximum sample available
    n_sample_win_max_orig = np.array([recovered_count, np.full_like(recovered_count, max_sample_size)]).min(axis=0)
    n_sample_win_max = n_sample_win_max_orig
    print("n_sample_win_max:", n_sample_win_max)
    print("n_remaining_sample:", n_remaining_sample)
    print("--------")

    while n_remaining_sample > 0:
        ## get the n_sample proportionally (this can exceed the recovered size)
        n_sample_win_tmp = np.random.choice(len(win_intervals), n_remaining_sample, p=recovered_freq)
        n_sample_win_tmp = np.bincount(n_sample_win_tmp, minlength=len(win_intervals))

        print("n_sample_win_tmp_org:", n_sample_win_tmp)
        print("n_sample_win_max:", n_sample_win_max)
        ## if sample_size is larger than the recovered size, make n_sample = n_sample_win_max
        n_sample_win_tmp = np.array([n_sample_win_tmp, n_sample_win_max]).min(axis=0)
        print("n_sample_win_tmp_after:", n_sample_win_tmp)


        ## update the number of available size remaining in the recovered individual & n_sample_win_final
        n_sample_win += n_sample_win_tmp
        recovered_count -= n_sample_win_tmp
        recovered_freq = recovered_count / recovered_count.sum()

        n_sample_win_max -= n_sample_win_tmp

        ## get the number of samples remaining
        n_remaining_sample = n_remaining_sample - n_sample_win_tmp.sum()


        print("recovered_count:", recovered_count)
        print("n_sample_win_max:", n_sample_win_max)
        print("n_remaining_sample:", n_remaining_sample)
        print("--------")

        if n_remaining_sample < 0:
            raise Exception('n_remaining_sample < 0')

    print("n_sample_win:", n_sample_win)
    print("n_sample_win_freq:", n_sample_win / n_sample_win.sum())

    # 3. sample from each window accordingly
    for this_win_idx in range(len(win_intervals)):
        print("this_win_idx:", this_win_idx, "/", len(win_intervals), n_sample_win[this_win_idx],
              recovered_count_orig[this_win_idx])
        this_win_end = win_intervals[this_win_idx]

        this_win_n_sample = n_sample_win[this_win_idx]

        if this_win_n_sample > 0:

            recovered_idx_win = np.where(recovered_window == this_win_end)[0]
            recovered_idx_long = np.repeat(recovered_idx_win, sim['recovered'][recovered_idx_win, 2].astype('int'))
            sampled_idx = np.random.choice(recovered_idx_long, size=this_win_n_sample)

            sampled_df = pd.DataFrame({'t': sim['recovered'][sampled_idx, 0],
                                       'genotype': sim['recovered'][sampled_idx, 1],
                                       'count': sim['recovered'][sampled_idx, 2],
                                       'clade': sim['recovered'][sampled_idx, 3]})

            clade_idx = 0
            clade_genotype_info = sim['genotype_info'][:, 2 * clade_idx: 2 * (clade_idx + 1)]

            if len(sampled_df) > 0:
                data['S_in_win'][this_win_idx, clade_idx] = statistics.n_segregating_from_samples(
                    sim['recovered'][sampled_idx, 1].astype(int),
                    # sampled_df['genotype'].to_numpy(dtype = 'int'),
                    clade_genotype_info)

                data['n_samples_in_win'][this_win_idx, clade_idx] = this_win_n_sample

    return data


def sampling (sim, params, data, dataset_dir, file_header, sampling_type,  n_sample=200, max_sampling_rate=None,max_sample_size = None, name_tail = "", imported_data = None):

    n_clade = sim['n_clade']
    state_var = sim['statevar']
    data['win_intervals'] = misc_.np_arange_prevent_imprecision(max(state_var[:, 0]), params['timeintro'], data['window_dt'] * (-1))[::-1]

    """
    if sampling_type == "_unif":
        sample_per_window = n_sample
        sample_name = dataset_dir + file_header + sampling_type + str(sample_per_window)
        data = sampledata.sampling_unif_fixed_with_SFS(sim, params, data, [sample_per_window], sample_name + name_tail + "_SFS.csv")
        
    if sampling_type == "_unif_maxrate":
        sample_per_window = n_sample
        sample_name = dataset_dir + file_header + sampling_type + str(sample_per_window) + "_" + str(
            max_sampling_rate).replace(".", "")
        data = sampledata.sampling_unif_fixed_maxrate_with_SFS(sim, params, data, sample_per_window, max_sampling_rate, sample_name + name_tail + "_SFS.csv")


    if sampling_type == "_full":
        sample_name = dataset_dir + file_header + sampling_type  # + "_" + str(n_sample_full)
        data, n_sample_full = sampledata.sampling_full_with_SFS(sim, params, data, sample_name + name_tail + "_SFS.csv")

   """


    if sampling_type == "_prop":
        n_total_sample = n_sample
        sample_name = dataset_dir + file_header + sampling_type + str(n_total_sample)
        data = sampling_prop_large_set(sim, data, n_total_sample)

    if sampling_type == "_data":
        imported_data = pd.DataFrame({'win_end': imported_data[:, 5],
                                      'n_sample': imported_data[:, 1],
                                      'clade': imported_data[:, 6]})

        n_total_sample = np.sum(imported_data['n_sample']).astype('int')
        sample_name = dataset_dir + file_header + sampling_type + str(n_total_sample)
        data = sampling_data_large_set(sim, data, imported_data)



    if sampling_type == "_prop_max":
        n_total_sample = n_sample
        sample_name = dataset_dir + file_header +  "_prop" +  str(n_total_sample) + "_max" + str(max_sample_size)
        data = sampling_prop_with_max(sim, data, n_total_sample, max_sample_size)



    np.savez(sample_name + name_tail+ ".npz", **{key: value for key, value in data.items()})
    dir_table = sample_name + name_tail+ ".tsv"
    f = open(dir_table, 'w')
    print('\t'.join(["window", "n", "s", "window_end", "clade"]), file=f)

    n_sample_sum = [0] * n_clade
    for i in range(len(data['win_intervals'])):
        for clade_idx in range(n_clade):
            data_values = [i, data['n_samples_in_win'][i][clade_idx], data['S_in_win'][i][clade_idx],  data['win_intervals'][i], clade_idx]
            print('\t'.join([str(x) for x in data_values]))
            print('\t'.join([str(x) for x in data_values]), file=f, flush=True)

            if ~np.isnan(data['n_samples_in_win'][i][clade_idx]):
                n_sample_sum[clade_idx] += data['n_samples_in_win'][i][clade_idx]


    print ("n_sample_sum: ", n_sample_sum)
    return





# ---------------------------
def _find_window(t_list, window):
    l = []
    for t in t_list:
        l.append(np.min(window[t <= window ], ))

    return np.array(l)

def _initiate_sample_data(data, win_intervals_size, n_clade):
    import copy
    sample_data = copy.deepcopy(data)

    sample_data['n_samples_in_win'] = np.full([win_intervals_size, n_clade], np.nan)
    sample_data['S_in_win'] = np.full([win_intervals_size, n_clade], np.nan)

    sample_data['sampled_individuals'] = np.empty([win_intervals_size, n_clade], dtype = 'object')


    return sample_data
