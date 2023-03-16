#!/bin/bash

PEM_PATH="$( cd "$( dirname "$0" )" && pwd -P )"
echo $PEM_PATH
# script below - - - - - -  - - - - - - -
france_data='gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv'
if [ $(pwd | grep ypar279 | wc -l) -eq 1 ];then   ## run from server
  case $1 in
  multi1224)
    simulname="france_multiple_te=191224"
    echo $1
    python $PEM_PATH/../eiuss/eiuss.py reconstruct \
          --config            empirical_data/france_data/configs/config_${simulname}.py \
          --input             empirical_data/france_data/data/${france_data} \
          --input_meanllk     empirical_data/france_data/${simulname}/gridsearch_log10eta,R0/meanllk.tsv \
          --outdir            empirical_data/france_data/${simulname} \
          --operators         log10eta R0 \
          --n_SMC_particles   200 \
          --n_grab            20 \
          --n_reps            1 \
          --n_save_particle   1 \
          --seed              2023${2}
    ;;


    multi0101)
      simulname="france_multiple_te=200101"
      echo $1
      python $PEM_PATH/../eiuss/eiuss.py reconstruct \
          --config            empirical_data/france_data/configs/config_${simulname}.py \
          --input             empirical_data/france_data/data/${france_data} \
          --input_meanllk     empirical_data/france_data/${simulname}/gridsearch_log10eta,R0/meanllk.tsv \
          --outdir            empirical_data/france_data/${simulname} \
          --operators         log10eta R0 \
          --n_SMC_particles   200 \
          --n_grab            20 \
          --n_reps            1 \
          --n_save_particle   1 \
          --seed              2023${2}
    ;;


    multi0108)
      simulname="france_multiple_te=200108"
      echo $1
      python $PEM_PATH/../eiuss/eiuss.py reconstruct \
          --config            empirical_data/france_data/configs/config_${simulname}.py \
          --input             empirical_data/france_data/data/${france_data} \
          --input_meanllk     empirical_data/france_data/${simulname}/gridsearch_log10eta,R0/meanllk.tsv \
          --outdir            empirical_data/france_data/${simulname} \
          --operators         log10eta R0 \
          --n_SMC_particles   200 \
          --n_grab            20 \
          --n_reps            1 \
          --n_save_particle   1 \
          --seed              2023${2}
    ;;

  esac

else      ## run from local
  case $1 in
  ready)
    sh ../run_remote.sh   bio_ready     empirical_data/france_data/run_reconstruct
    ;;

  run)
    sh ../run_remote.sh   bio_run       ../session_2023-03-13.log   'for i in {0..9}; do sbatch empirical_data/france_data/run_reconstruct_server.sh  multi1224 ${i}; done'
    sh ../run_remote.sh   bio_run       ../session_2023-03-13.log   'for i in {0..9}; do sbatch empirical_data/france_data/run_reconstruct_server.sh  multi0101 ${i}; done'
    sh ../run_remote.sh   bio_run       ../session_2023-03-13.log   'for i in {0..9}; do sbatch empirical_data/france_data/run_reconstruct_server.sh  multi0108 ${i}; done'
    ;;

  download)
    echo "-------------";
    echo "bio"
    sh ../run_remote.sh  bio_download   empirical_data/france_data/france_multiple_te=200108/reconstruct_log10eta,R0
    sh ../run_remote.sh  bio_download   empirical_data/france_data/france_multiple_te=200101/reconstruct_log10eta,R0
    sh ../run_remote.sh  bio_download   empirical_data/france_data/france_multiple_te=191224/reconstruct_log10eta,R0
    ;;


  test)
    for simulname in "france_multiple_te=191224" "france_multiple_te=200101" "france_multiple_te=200108";do
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
          --seed              202303

    done





  esac

fi
