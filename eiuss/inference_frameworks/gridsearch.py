import time
import importlib
import numpy as np
from itertools import product

import likelihood
# ------------------
def set_param_values(params, new_vals):
    for val_idx, val in enumerate(new_vals):
        if params['operators'][val_idx] in params.keys():
            params[params['operators'][val_idx]] = val
        else:
            raise Exception(
                'parameters to operate on must be defined in the params or params["model_params"] dictionary')
    return (params)


def run(params, imported_data, config):
    print(f'>>> epidemiological model used: {params["model_name"]}')
    particle_model_specific = importlib.import_module(params["model_name"] + "." + "particle")

    print(f'>>>> operating on the following parameters: {params["operators"]}')
    print(f'>>>> range for each parameter defined as: {list(zip(params["operators"], params["range"]))}')
    print(f'running {params["n_reps"]} replicates per parameter combination')

    ## make parameter combinations
    combos = list(product(*[np.arange(params['n_reps'])] + params['range']))

    ## make arrays to store parameter combinations and likelihood
    llh_out = np.zeros((len(combos), len(combos[0]) + 1))

    for combo_idx, combo in enumerate(combos):
        start = time.time()
        params = set_param_values(params, combo[1:])
        params = config.get_event_rate(params)

        ## initiate particles
        particle_init_size = np.ceil((imported_data['window_end'].max() - params['timestart']) / params['dt']).astype(
            int)
        print("particle_init_size: ", particle_init_size)
        particles = [particle_model_specific.Particle(1, params['n_class'],
                                                      init_size=particle_init_size + len(params['operators'])) for x in
                     range(params['n_SMC_particles'])]  # allocate particles

        ## get likelihood
        combo_llh, particles, _ = likelihood.calculate(particles, params, imported_data, config)
        print("combo_llh", combo_llh)

        llh_out[combo_idx, :len(combo)] = combo
        llh_out[combo_idx, -1] = combo_llh

        end = time.time()
        print(f'{end - start} seconds for parameter value/s {combo} and llh = {combo_llh}')

        ## save output
        import pickle
        print(params['out_name'])
        pickle.dump(llh_out, open(params['out_name'] + '.pkl', 'wb'))

    return (llh_out)