import os
import time
import pickle
import importlib
import numpy as np
import pandas as pd
from itertools import product

import utils
import likelihood


def change_param_values(params, new_vals):
    for val_idx, val in enumerate(new_vals):
        if params['operators'][val_idx] in params.keys():
            params[params['operators'][val_idx]] = val
        else:
            raise Exception('parameters to operate on must be defined in the params or params["model_params"] dictionary')
    return (params)


def run(params, imported_data, config):
    ## load input data
    print("\n" + "- " * 30)
    print(f'>>> segregating site trajectory data used: {params["input_data"]}', flush=True)
    imported_data, wind_dt = utils.import_data_simple(params["input_data"])
    params['window_dt'] = wind_dt
    print(f'>>>> window_dt inferred from the dataset : {params["window_dt"]}', flush=True)
    print(imported_data)

    ## load model-specific modules
    print("\n" + "- " * 30)
    print(f'>>> epidemiological model used: {params["model_name"]}', flush=True)
    #likelihood_model_specific = importlib.import_module("models." + params["model_name"] + "." + "likelihood")
    particle_model_specific   = importlib.import_module("models." + params["model_name"] + "." + "particle")


    ## make parameter combinations
    print("\n" + "- " * 30)
    print(f'>>>> operating on the following parameters: {params["operators"]}')
    print(f'>>>> range for each parameter defined as: {list(zip(params["operators"], params["range"]))}')
    print(f'     running {params["n_reps"]} replicates with random parameter combinations', flush=True)
    print("\n" + "- " * 30, flush=True)


    combos = list(product(*[np.arange(params['n_reps'])] + params['range']))


    if params['continued']:
        file_location = params['out_name'].replace(params['out_name'].rpartition("/")[-1], "")
        files_done = [x for x in os.listdir(file_location) \
                      if x.find(params['out_name'].rpartition("/")[-1]+ ".pkl") > -1 or
                        (x.find(params['out_name'].rpartition("/")[-1]+ "_continued") > -1 and x.find(".pkl") > -1)]
        # last_combo_idx = 0
        # for file in files_done:
        #     pkl = pickle.load(open(file_location + "/" + file, "rb"))
        #     last_combo_idx += np.where(pkl.sum(axis=1) == 0)[0].min()
        #
        # combos = combos[last_combo_idx:]

        last_combos = []
        for file in files_done:
            pkl = pickle.load(open(file_location + "/" + file, "rb"))
            last_combos += [x[:-1] for x in pkl.tolist() if not sum(x) == 0]

        combos = [list(x) for x in combos if not np.any([np.allclose(x,y) for y in last_combos])]
        combos = sorted(combos, key=lambda x: (x[2], x[1]), reverse=True)

        params['out_name'] += '_continued' + str(len([x for x in files_done if x.find('_continued') > -1]))
        print("\n" + "- " * 30)
        print(f'>>>> related files: {files_done}')
        print(f'>>>> continuing the gridsearch from combo: ')
        for i in combos:
            print(f'{i}')
        print("\n" + "- " * 30, flush=True)




    llh_out = np.zeros((len(combos), len(combos[0])+1))
    for combo_idx, combo in enumerate(combos):
        start = time.time()
        params = change_param_values(params, combo[1:])

        ## initiate particles
        particle_init_size = np.ceil((imported_data['window_end'].max()- params['timestart']) / params['dt']).astype(int)
        particles = [particle_model_specific.Particle(1, params['n_class'],  init_size= particle_init_size+ len(params['operators'])) for x in range(params['n_SMC_particles'])]  # allocate particles
        print(f"combo {combo_idx:7d}: {combo} (particle_init_size= {particle_init_size:7d})", flush=True)


        ## get likelihood
        combo_llh, particles, llh_window = likelihood.calculate(particles, params, imported_data, config)
        print("combo_llh", combo_llh)

        llh_out[combo_idx,:len(combo)]  = combo
        llh_out[combo_idx,-1]           = combo_llh

        end = time.time()
        print(f'{end - start} seconds for parameter value/s {combo} and llh = {combo_llh}', flush=True)


        # save output
        print(params['out_name'])
        pickle.dump(llh_out, open(params['out_name'] +'.pkl', 'wb'))
        pd.DataFrame(llh_out, columns=["n_iter"] + params["operators"] + ["llh"]).\
            to_csv(params['out_name'] + '_loglikelihood'  + '.tsv', sep = "\t", index=False)

    return(llh_out)