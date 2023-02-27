# ------------------
## Epidemiological Inference using Segregating Sites
## last modified : 2023-01-24 by Yeongseon Park
## contact: yeongseon.park@emory.edu
# ------------------
import sys
import os
import time
import pickle
import shutil
import argparse
import numpy as np
from importlib import import_module

import tool


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config')
    parser.add_argument('--inputData')
    args = parser.parse_args()

    if args.inputData == None:
        print('>>>> no input data provided, using example dataset')
        #args.inputData = "../data_simulated/dataset_piuss_seed1234/SEIR_simple_prop150_sampled_during_28-56_1219.tsv"
        #args.inputData = "../data_France/empty.tsv"
        args.inputData = "../data_France/mm_0518/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_match_largest_s.tsv"


    if args.config == None:
        print('>>>> no input config file provided, using example config file (SEIR)')
        #args.config = "config/config_sim-early45_2Dprofilelikelihood.py"
        #args.config = "config/config_hmSE1E2IR_simulate.py"
        #args.config = "model_hmSE1E2IR/config_hmSE1E2IR_simulate.py"
        #args.config = "config/config_France_SE1E2IR_hetero_migration.py"
        args.config = "config/config_France_multiple_exp_te=191224_2Dgridsearch.py"


    # load in input data
    print(f'>>>> importing dataset from {args.inputData}')
    imported_data = tool.import_data_simple(args.inputData)
    print(imported_data)

    # load in config file
    print (f'>>>> import model configureation from {args.config}', flush=True)
    config = import_module(args.config.replace('.py', '').replace('/', '.'))
    params = config.set_param(imported_data)

    os.makedirs(params['out_name'].rpartition("/")[0], exist_ok=True)
    shutil.copy(args.config, params['out_name'].rpartition("/")[0])




    # now, send off to inference framework
    inference_module_name = params['inference_framework'].replace('.py', '').replace('/', '.')
    print(f'>>>> using inference framwork: {inference_module_name}')
    inference_module_name = 'inference_frameworks.' + inference_module_name
    inference = import_module(inference_module_name)

    out_dat = inference.run(params, imported_data, config)
    print(params['out_name'])
    pickle.dump(out_dat, open(params['out_name']+'.pkl', 'wb'))
    
    

    print(f'>>>> inference complete')


if __name__ == "__main__":
    run()