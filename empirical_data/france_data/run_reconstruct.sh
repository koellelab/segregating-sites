#!/bin/bash

PEM_PATH="$( cd "$( dirname "$0" )" && pwd -P )"
echo $PEM_PATH
# script below - - - - - -  - - - - - - -
france_data='gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv'
case $1 in
  multi1224)
    echo $1
    simulname="france_multiple_te=191224"
    python $PEM_PATH/../../eiuss/eiuss.py reconstruct \
        --config            empirical_data/france_data/configs/config_${simulname}.py \
        --input             empirical_data/france_data/data/${france_data} \
        --input_meanllk     empirical_data/france_data/${simulname}/gridsearch_log10eta,R0/meanllk.tsv \
        --outdir            empirical_data/france_data/${simulname} \
        --operators         log10eta R0 \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            10 \
        --n_save_particle   1 \
        --seed              202302
    ;;


    multi0101)
    echo $1
    simulname="france_multiple_te=200101"
    python $PEM_PATH/../../eiuss/eiuss.py reconstruct \
        --config            empirical_data/france_data/configs/config_${simulname}.py \
        --input             empirical_data/france_data/data/${france_data} \
        --input_meanllk     empirical_data/france_data/${simulname}/gridsearch_log10eta,R0/meanllk.tsv \
        --outdir            empirical_data/france_data/${simulname} \
        --operators         log10eta R0 \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            10 \
        --n_save_particle   1 \
        --seed              202302
    ;;


    multi0108)
    echo $1
    simulname="france_multiple_te=200108"
    python $PEM_PATH/../eiuss/eiuss.py reconstruct \
        --config            empirical_data/france_data/configs/config_${simulname}.py \
        --input             empirical_data/france_data/data/${france_data} \
        --input_meanllk     empirical_data/france_data/${simulname}/gridsearch_log10eta,R0/meanllk.tsv \
        --outdir            empirical_data/france_data/${simulname} \
        --operators         log10eta R0 \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            10 \
        --n_save_particle   1 \
        --seed              202302
    ;;

    esac