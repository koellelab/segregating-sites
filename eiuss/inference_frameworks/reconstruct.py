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
    print(f'>>> segregating site trajectory data used: {params["input_data"]}')
    imported_data, wind_dt = utils.import_data_simple(params["input_data"])
    params['window_dt'] = wind_dt
    print(f'>>>> window_dt inferred from the dataset : {params["window_dt"]}')
    print(imported_data)

    ## load model-specific modules
    print("\n" + "- " * 30)
    print(f'>>> epidemiological model used: {params["model_name"]}')
    #likelihood_model_specific = importlib.import_module("models." + params["model_name"] + "." + "likelihood")
    particle_model_specific   = importlib.import_module("models." + params["model_name"] + "." + "particle")


    ## make parameter combinations
    print("\n" + "- " * 30)
    print(f'>>>> parameter combinations with top 10 llk : {params["input_meanllk"]}')
    df_mean_llk = pd.read_csv(params['input_meanllk'], sep="\t")
    llk = np.exp(df_mean_llk.mean_logL) / np.exp(df_mean_llk.mean_logL).sum()
    llk = np.where(~np.isnan(llk), llk, 0)
    print(df_mean_llk.sort_values('mean_logL', ascending=False).iloc[:10].to_string())

    ## choose top 10
    df_mean_llk = df_mean_llk.sort_values('mean_logL', ascending=False)
    sampled_idx = np.arange(params["n_reps"])
    #sampled_idx = params['rng'].choice(len(df_mean_llk.mean_logL), size=params["n_reps"], p= llk/llk.sum(), replace=True)
    sampled_parameter = df_mean_llk.iloc[sampled_idx][params["operators"]]
    combos = [tuple([idx] + sampled_parameter.iloc[idx].values.tolist()) for idx in range(len(sampled_parameter))]

    print("\n" + "- " * 30)
    print(f'>>>> reconstruct state variables with from combo:{combos}')
    print(df_mean_llk.iloc[sampled_idx].to_string())
    print("\n" + "- " * 30)


    llh_out = np.zeros((len(combos), len(combos[0])+1))
    for combo_idx, combo in enumerate(combos):
        start = time.time()
        params = change_param_values(params, combo[1:])

        ## initiate particles
        particle_init_size = np.ceil((imported_data['window_end'].max()- params['timestart']) / params['dt']).astype(int)
        particles = [particle_model_specific.Particle(1, params['n_class'],  init_size= particle_init_size+ len(params['operators'])) for x in range(params['n_SMC_particles'])]  # allocate particles
        print(f"combo {combo_idx:7d}: {combo} (particle_init_size= {particle_init_size:7d})")


        ## get likelihood
        combo_llh, particles, llh_window = likelihood.calculate(particles, params, imported_data, config)
        print("combo_llh", combo_llh)

        llh_out[combo_idx,:len(combo)]  = combo
        llh_out[combo_idx,-1]           = combo_llh

        end = time.time()
        print(f'{end - start} seconds for parameter value/s {combo} and llh = {combo_llh}', flush=True)

        # save output
        print(params['out_name'])
        pickle.dump(llh_out, open(params['out_name'] + '.pkl', 'wb'))
        pd.DataFrame(llh_out, columns=["n_iter"] + params["operators"] + ["llh"]). \
            to_csv(params['out_name'] + '_loglikelihood' + '.tsv', sep="\t", index=False)

        # save particle
        save_idx = params['rng'].choice(len(particles), size=params['n_save_particle'])

        if hasattr(particles[0], "infected_from_outside"):
            saved_particles = [[particles[i].statevar, particles[i].n_segregating, particles[i].n_genotype, particles[i].n_site] for i in save_idx]
        else:
            saved_particles = [
                [particles[i].statevar, particles[i].n_segregating, particles[i].n_genotype, particles[i].n_site, particles[i].infected_from_outside] for i in save_idx]

        pickle.dump(saved_particles, open(params['out_name'] + str(combo).replace(")", f", {df_mean_llk.iloc[sampled_idx[combo_idx]]['mean_logL']:.3f})") +  '.pkl', 'wb'))




    return (llh_out)