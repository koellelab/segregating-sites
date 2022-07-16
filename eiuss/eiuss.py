# ------------------------
# < Epidemiological inference for
#  emerging viruses using segregating sites >
#  Yeongseon Park, Michael Martin, Katia Koelle
#
# ------------------------

import os
import pickle
import shutil
import argparse
import numpy as np
import pandas as pd
from importlib import import_module


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config')
    parser.add_argument('--inputData')
    args = parser.parse_args()

    if args.inputData == None:
        print('>>>> no input data provided, using example dataset')
        args.inputData = "../data_example/example_data.tsv"

    if args.config == None:
        print('>>>> no input config file provided, using example config file (SEIR)')
        args.config = "config/config_example.py"

    # load in input data
    print(f'>>>> importing dataset from {args.inputData}')
    imported_data = pd.read_csv(args.inputData, sep='\t')
    imported_data = imported_data[imported_data['n'].notnull()]  # trim nan rows

    # load in config file
    print(f'>>>> import model configureation from {args.config}')
    config = import_module(args.config.replace('.py', '').replace('/', '.'))
    params = config.set_param(imported_data)

    # make folder to store output files
    out_folder = params['out_name'].rpartition("/")[0]
    os.makedirs(out_folder, exist_ok=True)
    shutil.copy(args.config, out_folder)

    # now, send off to inference framework
    inference_module_name = params['inference_framework'].replace('.py', '').replace('/', '.')
    print(f'>>>> using inference framwork: {inference_module_name}')
    inference_module_name = 'inference_frameworks.' + inference_module_name
    inference = import_module(inference_module_name)

    # run analysis and save output
    out_dat = inference.run(params, imported_data, config)
    print(params['out_name'])
    pickle.dump(out_dat, open(params['out_name'] + '.pkl', 'wb'))

    print(f'>>>> inference complete')


if __name__ == "__main__":
    run()