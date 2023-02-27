import os
import sys
import pickle
import shutil
import argparse
import numpy as np
from datetime import date
from importlib import import_module

import pandas as pd

eiuss_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(eiuss_dir)
# ---------------------------
from datetime import datetime
from datetime import timedelta

def convert_to_matlab_date (datestr, sep = "/"):
    return datetime.strptime(datestr, sep.join(['%Y', '%m', '%d'])).toordinal() + 366

def convert_to_caldate (matlab_datenum, sep = "/"):
    return (datetime(1, 1, 1) + timedelta(matlab_datenum - 367)).strftime(sep.join(['%Y', '\n%m', '%d']))



def convert_range_gridsearch(str_range):
    if str_range.find(":"):
        parsed = [convert_to_matlab_date(x) if x.find("/") > -1 else float(x) for x in str_range.split(":")]
        return np.arange(parsed[0], parsed[1], parsed[2])

    elif str_range.find(","):
        return np.array([float(x) for x in str_range.split(",")])


def convert_range_randsearch(str_range):
    if str_range.find(":"):
        parsed = [float(x) for x in str_range.split(":")]
        return [parsed[0], parsed[1]]


def run():
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(dest='inference_methods')

    simul_epi  = subparser.add_parser('simulate_epidemic')
    simul_samp = subparser.add_parser('simulate_sampling')
    gridsearch = subparser.add_parser('gridsearch')
    randomsearch = subparser.add_parser('randomsearch')
    reconstruct = subparser.add_parser('reconstruct')
    log_particlefiltering = subparser.add_parser('log_particlefiltering')

    simul_epi.add_argument('--config', type=str, required=True)
    simul_epi.add_argument('--type', type=str, required=True)
    simul_epi.add_argument('--seed',   type=int, required=False,    default=int(date.today().strftime("%Y%m%d")))
    simul_epi.add_argument('--n_trial', type=int,                   default=1)

    simul_samp.add_argument('--config', type=str, required=True)
    simul_samp.add_argument('--input', type=str, required=True)
    simul_samp.add_argument('--input_samp', type=str, required=True)
    simul_samp.add_argument('--sample_type', type=str, required=True)
    simul_samp.add_argument('--sample_size', type=int)
    simul_samp.add_argument('--window_dt', type=int, required=True)
    simul_samp.add_argument('--sample_between', type=int, nargs=2,  default=None)
    simul_samp.add_argument('--seed', type=int, required=False,     default=int(date.today().strftime("%Y%m%d")))
    simul_samp.add_argument('--n_trial', type=int,                  default=1)

    gridsearch.add_argument('--config', type=str, required=True)
    gridsearch.add_argument('--input', type=str, required=True)
    gridsearch.add_argument('--seed', type=int, required=False,     default=int(date.today().strftime("%Y%m%d")))
    gridsearch.add_argument('--n_SMC_particles', type=int)
    gridsearch.add_argument('--n_grab', type=int)
    gridsearch.add_argument('--n_reps', type=int)
    gridsearch.add_argument('--operators', type=str, nargs='+')
    gridsearch.add_argument('--range', type=str)
    gridsearch.add_argument('--n_trial', type=int, default=1)
    gridsearch.add_argument('--outdir', type=str, default=None)
    gridsearch.add_argument('--continued', action='store_true')

    randomsearch.add_argument('--config', type=str, required=True)
    randomsearch.add_argument('--input', type=str, required=True)
    randomsearch.add_argument('--seed', type=int, required=False,     default=int(date.today().strftime("%Y%m%d")))
    randomsearch.add_argument('--n_SMC_particles', type=int)
    randomsearch.add_argument('--n_grab', type=int)
    randomsearch.add_argument('--n_reps', type=int)
    randomsearch.add_argument('--operators', type=str, nargs='+')
    randomsearch.add_argument('--range', type=str)
    randomsearch.add_argument('--n_trial', type=int, default=1)


    reconstruct.add_argument('--config', type=str, required=True)
    reconstruct.add_argument('--input', type=str, required=True)
    reconstruct.add_argument('--input_meanllk', type=str, required=True)
    reconstruct.add_argument('--seed', type=int, required=False,     default=int(date.today().strftime("%Y%m%d")))
    reconstruct.add_argument('--operators', type=str, nargs='+')
    reconstruct.add_argument('--n_SMC_particles', type=int)
    reconstruct.add_argument('--n_grab', type=int)
    reconstruct.add_argument('--n_reps', type=int)
    reconstruct.add_argument('--n_save_particle', type=int)
    reconstruct.add_argument('--n_trial', type=int, default=1)


    log_particlefiltering.add_argument('--config', type=str, required=True)
    log_particlefiltering.add_argument('--input', type=str, required=True)
    log_particlefiltering.add_argument('--seed', type=int, required=False,     default=int(date.today().strftime("%Y%m%d")))
    log_particlefiltering.add_argument('--n_SMC_particles', type=int)
    log_particlefiltering.add_argument('--n_grab', type=int)
    log_particlefiltering.add_argument('--operators', type=str, nargs='+')
    log_particlefiltering.add_argument('--values', type=float, nargs='+')
    log_particlefiltering.add_argument('--outdir', type=str, default=None)
    log_particlefiltering.add_argument('--continued', action='store_true')
    log_particlefiltering.add_argument('--n_trial', type=int, default=1)
    log_particlefiltering.add_argument('--n_reps', type=int, default=1)

    args = parser.parse_args()

    # load in input data


    ## load in config file
    shutil.copy(args.config, f'{eiuss_dir}/config/{args.config.rpartition("/")[-1]}')
    print (f'>>>> import model configureation from {args.config}', flush=True)
    config = import_module(f'config.{args.config.rpartition("/")[-1]}'.replace(".py", ""))
    params = config.set_param(args.seed)


    ## make output dir
    if args.inference_methods == 'simulate_epidemic':
        output_dir = args.config.replace("configs/", "").replace("config_", "").replace(".py", "")

        params['out_name'] = f'{output_dir}/seed{args.seed}'
        imported_data = None
        args.inference_methods += "s_" + args.type

    if args.inference_methods == 'simulate_sampling':
        if args.sample_type == "full":
            args.sample_size = ''
            args.sample_between = None

        elif args.sample_type == "reset_window":
            args.sample_size = None
            params['input_samp'] = args.input_samp


        if args.sample_type == "reset_window":
            output_dir = args.input_samp.replace(args.input_samp.split("/")[-2],args.input_samp.split("/")[-2] + "_resetwin" + str(args.window_dt))
        else:
            output_dir = f'{args.input}_{args.sample_type}{args.sample_size}_win{args.window_dt}'

        if args.sample_between:
            output_dir += "_" + "-".join([str(x) for x in args.sample_between])
        print(">> ", output_dir)

        params = vars(args)
        params['out_name'] = f'{output_dir}/seed{args.seed}'
        params['rng']      = np.random.default_rng(args.seed)

        config        = None
        imported_data = None        ## currently, not supported

    if args.inference_methods == 'gridsearch':
        print(args.range)
        args.range = [convert_range_gridsearch(x) for x in args.range.strip().split(";")]

        if args.outdir:
            output_dir = args.outdir + f'/gridsearch_{",".join(args.operators)}'
        else:
            output_dir = f'{args.input.replace("_segsites.tsv", "")}_gridsearch_{",".join(args.operators)}'

        print(output_dir)
        params['out_name']   = f'{output_dir}/seed{args.seed}'
        params["input_data"] = args.input
        params["operators"] = args.operators
        params["range"] = args.range
        params["n_grabs"] = args.n_grab
        params["n_reps"] = args.n_reps
        params["n_SMC_particles"] = args.n_SMC_particles
        params["continued"] = args.continued


        imported_data = None  ## currently, not supported

    if args.inference_methods == 'randomsearch':
        args.range = [convert_range_randsearch(x) for x in args.range.split(";")]
        print(args.range)
        output_dir = f'{args.input.replace("_segsites.tsv", "")}_randomsearch_{",".join(args.operators)}'

        print(">>seed" ,  args.seed, flush=True)
        params['out_name']   = f'{output_dir}/seed{args.seed}'
        params["input_data"] = args.input
        params["operators"] = args.operators
        params["range"] = args.range
        params["n_grabs"] = args.n_grab
        params["n_reps"] = args.n_reps
        params["n_SMC_particles"] = args.n_SMC_particles

        imported_data = None  ## currently, not supported

    if args.inference_methods == 'reconstruct':
        output_dir = f'{args.input.replace("_segsites.tsv", "")}_reconstruct'

        print(">>seed", args.seed, flush=True)
        params['out_name']   = f'{output_dir}/seed{args.seed}'
        params["input_data"] = args.input
        params["input_meanllk"] = args.input_meanllk
        params["operators"] = args.operators
        params["n_SMC_particles"] = args.n_SMC_particles
        params["n_grabs"] = args.n_grab
        params["n_reps"] = args.n_reps
        params["n_save_particle"] = args.n_save_particle

        imported_data = None  ## currently, not supported





    if args.inference_methods == 'log_particlefiltering':
        if args.outdir:
            output_dir = args.outdir + f'/showparticles_{",".join(args.operators)}'
        else:
            output_dir = f'{args.input.replace("_segsites.tsv", "")}_showparticles_{",".join(args.operators)}'

        print(output_dir)
        params['out_name']   = f'{output_dir}/seed{args.seed}'
        params["input_data"] = args.input
        params["operators"] = args.operators
        params["range"] = [[x] for x in args.values]
        params["n_reps"] = args.n_reps
        params["n_grabs"] = args.n_grab
        params["n_SMC_particles"] = args.n_SMC_particles
        params["continued"] = args.continued

        imported_data = None  ## currently, not supported



    os.makedirs(output_dir, exist_ok=True)

    # now, send off to inference framework
    params['out_name_orig'] = params['out_name']
    for idx_trial in range(args.n_trial):
        if args.n_trial == 1:
            params['out_name'] = params['out_name_orig']
        else:
            params['out_name'] = params['out_name_orig']+f'_trial{idx_trial}'
            #params['out_name'] = params['out_name_orig'] + f'_trial{idx_trial:0{int(np.log10(args.n_trial)) + 1}}'

        inference_module_name = args.inference_methods
        print(f'>>>> using inference framwork: {inference_module_name}')
        inference_module_name = 'inference_frameworks.' + inference_module_name
        inference = import_module(inference_module_name)

        out_dat = inference.run(params, imported_data, config)
        print(params['out_name'])
        #pickle.dump(out_dat, open(params['out_name']+'.pkl', 'wb'))

        #params['seed'] = args.seed + 1
        #params['rng']  = np.random.default_rng(args.seed + 1)
        #params['rng2'] = np.random.default_rng(args.seed + 1)


    os.remove(f'{eiuss_dir}/config/{args.config.rpartition("/")[-1]}')

    return 


if __name__ == "__main__":
    run()