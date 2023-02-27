#!/bin/bash
########################################
# 2023-01-27
# script for Figure 3 and 5
# calculate the likelihood for a range of (R0, t0) combinations
# : rng for mutations and epidemics are separated
#   - by setting param['rng2'] = param['rng'],
#     it works like an older version
########################################
PEM_PATH="$( cd "$( dirname "$0" )" && pwd -P )"
echo $PEM_PATH
# script below - - - - - -  - - - - - - -

simulname="simpleSEIR_R0=16e-1_mu2e-1"
dataname="seed1234_prop500_win4/seed230201"

case $1 in
  reconstruct)
    echo $1
    python $PEM_PATH/../eiuss/eiuss.py reconstruct \
        --config            simulated_data/configs/config_${simulname}.py \
        --input             simulated_data/${simulname}/${dataname}_segsites.tsv \
        --input_meanllk     simulated_data/${simulname}/${dataname}_gridsearch_R0,timestart/meanllk.tsv \
        --operators         R0 timestart\
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            10 \
        --n_save_particle   1 \
        --seed              202302${3}\
    ;;

  particle)
    echo $1
    python $PEM_PATH/../eiuss/eiuss.py log_particlefiltering \
        --config            simulated_data/configs/config_${simulname}.py \
        --input             simulated_data/${simulname}/${dataname}_segsites.tsv \
        --operators         R0 timestart\
        --values            1.7 16 \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --seed              202302${3}\

    ;;


  esac

