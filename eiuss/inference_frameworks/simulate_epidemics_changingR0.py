import sys
import copy
import time
import importlib
import numpy as np
import pandas as pd

import utils

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# ------------------
def set_param_values(params, new_vals):
    for val_idx, val in enumerate(new_vals):
        if params['operators'][val_idx] in params.keys():
            params[params['operators'][val_idx]] = val
        else:
            raise Exception('parameters to operate on must be defined in the params or params["model_params"] dictionary')
    return (params)





def run (params, imported_data, config):

    print(f'>>> epidemiological model used: {params["model_name"]}')
    particle_model_specific = importlib.import_module("models." + params["model_name"] + "." + "particle")

    sim_window = utils.np_arange_prevent_imprecision(params['timeend'], params['timestart'],-1 * params['dt'])[::-1]
    while True:
        time_start = time.time()
        particles = [particle_model_specific.Particle(1, params['n_class'], len(sim_window)+1)]
        particle  = particles[0]

        # initialize particles & add n_segregating row if needed
        particle.initiate(params, len(sim_window))
        particle.n_segregating[:len(sim_window), -1] = sim_window
        particle.recovered = pd.DataFrame()
        not_yet_updated = True


        this_time_start = params['timestart']
        for idx in range(len(sim_window)):
            this_time_end = sim_window[idx]
            # simulate one step
            params['model_nextdt'](particle, params, this_time_start, this_time_end)
            # record the genotype recovered during this step
            sample_genotype     = particle.genotype_count[:, params['idx_sampling']] > 0
            particle.recovered  = particle.recovered.append(
                pd.DataFrame({'t'       : this_time_end,
                              'genotype': np.where(sample_genotype)[0].astype('int'),
                              'count'   : particle.genotype_count[:, params['idx_sampling']][sample_genotype].astype('int')}))


            ## change R0 when reaching the threshold for the first time
            if (particle.statevar[particle.last_idx, params['idx_virus']+2].sum() >  params['threshold']) and not_yet_updated:
                params['R0'] = params['reduced_R0']
                params = config.update_dependent_param(params)
                not_yet_updated = False


            this_time_start = this_time_end

        printStr = f"n_genotype: {particle.n_genotype}, n_site: {particle.n_site}\n"
        print(printStr)

        if particle.n_genotype > 100:
            print("runtime: ", time.time() - time_start, "(sec)")
            particle.trim_arrays()
            particle.recovered = particle.recovered.to_numpy()

            np.savez  (params['out_name'] + "_simtaj.npz", **{key: np.array(value) for key, value in particle.__dict__.items() if key != 'model_nextdt'})
            np.savez  (params['out_name'] + "_params.npz", **{key: np.array(value) for key, value in params.items() if key != 'model_nextdt'})
            np.savetxt(params['out_name'] + '_statevar.tsv', particle.statevar, delimiter = "\t", header = '\t'.join(["t", "S"] + params['compartments']))

            break


    return



















