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
if [ $(pwd | grep ypar279 | wc -l) -eq 1 ];then   ## run from server
  case $1 in
  full)
    echo full
    ## full dataset with mu = 0.2
    simulname="simpleSEIR_R0=16e-1_mu2e-1"
    dataname="seed1234_prop500_win4/seed230201"
    python $PEM_PATH/eiuss/eiuss.py gridsearch \
        --config            simulated_data/configs/config_${simulname}.py \
        --input             simulated_data/${simulname}/${dataname}_segsites.tsv \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            1 \
        --operators         R0 timestart\
        --range             "1.${2}:2.51:0.2;30:-50.1:-2"\
        --seed              202302${3}\
        ${4}
      ;;

  early_low)
    echo early_low
    ## early dataset with mu = 0.2
    simulname="simpleSEIR_R0=16e-1_mu2e-1"
    dataname="seed1234_unif10_win4_32-52/seed230201"
    python $PEM_PATH/eiuss/eiuss.py gridsearch \
        --config            simulated_data/configs/config_${simulname}.py \
        --input             simulated_data/${simulname}/${dataname}_segsites.tsv \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            1 \
        --operators         R0 timestart\
        --range             "1.${2}:6.01:0.2;36:-50.1:-2"\
        --seed              202302${3}\
        ${4}
      ;;

  early_high)
    echo early_high
    ## early dataset with mu = 0.2
    simulname="simpleSEIR_R0=16e-1_mu4e-1"
    dataname="seed221110_unif10_win4_32-52/seed230201"
    python $PEM_PATH/eiuss/eiuss.py gridsearch \
        --config            simulated_data/configs/config_${simulname}.py \
        --input             simulated_data/${simulname}/${dataname}_segsites.tsv \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            1 \
        --operators         R0 timestart\
        --range             "1.${2}:6.01:0.2;36:-50.1:-2"\
        --seed              202302${3}\
        ${4}
      ;;

  full_extend)
    echo full
    ## full dataset with mu = 0.2
    simulname="simpleSEIR_R0=16e-1_mu2e-1"
    dataname="seed1234_prop500_win4/seed230201"
    python $PEM_PATH/eiuss/eiuss.py gridsearch \
        --config            simulated_data/configs/config_${simulname}.py \
        --input             simulated_data/${simulname}/${dataname}_segsites.tsv \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            1 \
        --operators         R0 timestart\
        --range             "1.${2}:2.51:0.2;36:30.1:-2"\
        --seed              202302${3}\
        ${4}
      ;;

  early_low_extend)
    echo early_low_extend
    ## early dataset with mu = 0.2
    simulname="simpleSEIR_R0=16e-1_mu2e-1"
    dataname="seed1234_unif10_win4_32-52/seed230201"
    python $PEM_PATH/eiuss/eiuss.py gridsearch \
        --config            simulated_data/configs/config_${simulname}.py \
        --input             simulated_data/${simulname}/${dataname}_segsites.tsv \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            20 \
        --operators         R0 timestart\
        --range             "1.0:6.01:0.1;${2}:${2}.1:-2"\
        --seed              202302${3}\
        ${4}
      ;;

  early_high_extend)
    echo early_high_extend
    ## early dataset with mu = 0.2
    simulname="simpleSEIR_R0=16e-1_mu4e-1"
    dataname="seed221110_unif10_win4_32-52/seed230201"
    python $PEM_PATH/eiuss/eiuss.py gridsearch \
        --config            simulated_data/configs/config_${simulname}.py \
        --input             simulated_data/${simulname}/${dataname}_segsites.tsv \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            20 \
        --operators         R0 timestart\
        --range             "1.0:6.01:0.1;${2}:${2}.1:-2"\
        --seed              202302${3}\
        ${4}
      ;;
  esac

else      ## run from local

  ## run analysis on 2023-Feb-02
  ## total 20 reps -> even numbers: 10 jobs in  bio cluster
  ##               ->             : 10 jobs in rsph cluster
  ##               ->  odd numbers: 20 jobs in rsph cluster

  case $1 in
  ready)
    sh run_remote.sh   bio_ready simulated_data/run_2D-loglikelihood
    #sh run_remote.sh  rsph_ready simulated_data/run_2D-loglikelihood
    ;;

  run)
    #sh run_remote.sh   bio_run session_2023-02-12.log 'for i in {0..9}    ;do sbatch simulated_data/run_2D-loglikelihood_server.sh full       0 ${i}; done'
    #sh run_remote.sh  rsph_run session_2023-02-12.log 'for i in {10..19}  ;do sbatch simulated_data/run_2D-loglikelihood_server.sh full       0 ${i}; done'
    #sh run_remote.sh  rsph_run session_2023-02-12.log 'for i in {20..29}  ;do sbatch simulated_data/run_2D-loglikelihood_server.sh full       1 ${i}; done'
    #sh run_remote.sh  rsph_run session_2023-02-12.log 'for i in {30..39}  ;do sbatch simulated_data/run_2D-loglikelihood_server.sh full       1 ${i}; done'
    #sh run_remote.sh  rsph_run session_2023-02-12.log 'for i in {40..49}  ;do sbatch simulated_data/run_2D-loglikelihood_server.sh full_extend 0 ${i}; done'    ## extend t0 up to 34
    #sh run_remote.sh  rsph_run session_2023-02-12.log 'for i in {50..59}  ;do sbatch simulated_data/run_2D-loglikelihood_server.sh full_extend 1 ${i}; done'    ## extend t0 up to 34

    #sh run_remote.sh  bio_run session_2023-02-12.log 'for i in {40..49}  ;do sbatch simulated_data/run_2D-loglikelihood_server.sh full_extend 0 ${i} --continued; done'    ## extend t0 up to 34
    #sh run_remote.sh  bio_run session_2023-02-12.log 'for i in {50..59}  ;do sbatch simulated_data/run_2D-loglikelihood_server.sh full_extend 1 ${i} --continued; done'    ## extend t0 up to 34
    #sh run_remote.sh  bio_run session_2023-02-12.log 'for i in {60..69}  ;do sbatch simulated_data/run_2D-loglikelihood_server.sh full_extend 0 ${i}; done'    ## extend t0 up to 34
    #sh run_remote.sh  bio_run session_2023-02-12.log 'for i in {70..79}  ;do sbatch simulated_data/run_2D-loglikelihood_server.sh full_extend 1 ${i}; done'    ## extend t0 up to 34


    sh run_remote.sh   bio_run  session_2023-02-12.log 'for i in {0..9}; do sbatch simulated_data/run_2D-loglikelihood_server.sh  early_low  0 ${i} --continued; done'
    sh run_remote.sh   bio_run  session_2023-02-12.log 'for i in {0..9}; do sbatch simulated_data/run_2D-loglikelihood_server.sh  early_high 0 ${i} --continued; done'
    sh run_remote.sh   bio_run  session_2023-02-12.log 'for i in {10..19}; do sbatch simulated_data/run_2D-loglikelihood_server.sh  early_low  0 ${i} --continued; done'
    sh run_remote.sh   bio_run  session_2023-02-12.log 'for i in {10..19}; do sbatch simulated_data/run_2D-loglikelihood_server.sh  early_high 0 ${i} --continued; done'
    sh run_remote.sh   bio_run  session_2023-02-12.log 'for i in {20..29} ;do sbatch simulated_data/run_2D-loglikelihood_server.sh  early_low  1 ${i} --continued; done'
    sh run_remote.sh   bio_run  session_2023-02-12.log 'for i in {20..29} ;do sbatch simulated_data/run_2D-loglikelihood_server.sh  early_high 1 ${i} --continued; done'
    sh run_remote.sh   bio_run  session_2023-02-12.log 'for i in {30..39} ;do sbatch simulated_data/run_2D-loglikelihood_server.sh  early_low  1 ${i} --continued; done'
    sh run_remote.sh   bio_run  session_2023-02-12.log 'for i in {30..39} ;do sbatch simulated_data/run_2D-loglikelihood_server.sh  early_high 1 ${i} --continued; done'

    #sh run_remote.sh   bio_run session_2023-02-17.log 'for i in {52..60..2};do sbatch simulated_data/run_2D-loglikelihood_server.sh   early_low_extend  -${i} ${i}; done'    ## extend t0 down to -60
    #sh run_remote.sh   bio_run session_2023-02-17.log 'for i in {52..60..2};do sbatch simulated_data/run_2D-loglikelihood_server.sh   early_low_extend  -${i} ${i}; done'    ## extend t0 down to -60
    #sh run_remote.sh   bio_run session_2023-02-17.log 'for i in {52..60..2};do sbatch simulated_data/run_2D-loglikelihood_server.sh   early_high_extend  -${i} ${i}; done'    ## extend t0 down to -60
    #sh run_remote.sh   bio_run session_2023-02-17.log 'for i in {52..60..2};do sbatch simulated_data/run_2D-loglikelihood_server.sh   early_high_extend  -${i} ${i}; done'    ## extend t0 down to -60
    ;;


  download)
    sh run_remote.sh  bio_download   simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_gridsearch_R0,timestart
    sh run_remote.sh  bio_download   simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/seed230201_gridsearch_R0,timestart
    sh run_remote.sh  bio_download   simulated_data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_unif10_win4_32-52/seed230201_gridsearch_R0,timestart

    sh run_remote.sh rsph_download   simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_gridsearch_R0,timestart
    sh run_remote.sh rsph_download   simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/seed230201_gridsearch_R0,timestart
    sh run_remote.sh rsph_download   simulated_data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_unif10_win4_32-52/seed230201_gridsearch_R0,timestart
    ;;

  esac
fi


