#!/bin/bash

PEM_PATH="$( cd "$( dirname "$0" )" && pwd -P )"
echo $PEM_PATH
# script below - - - - - -  - - - - - - -
if [ $(pwd | grep ypar279 | wc -l) -eq 1 ];then   ## run from server
  case $1 in
  single)
    echo $1 "| seed = " ${4}
    ## full dataset with mu = 0.2
    simulname="france_single"
    python $PEM_PATH/../eiuss/eiuss.py gridsearch \
        --config            empirical_data/france_data/configs/config_${simulname}.py \
        --input             empirical_data/france_data/data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
        --outdir            empirical_data/france_data/${simulname} \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            1 \
        --operators         timestart R0 \
        --range             "2020/02/17:2019/11/30:-3;${2}:${3}:0.1" \
        --seed              202303${4}\
        ${5}
      ;;

  multi1224)
    echo $1 "| seed = " ${4}
    simulname="france_multiple_te=191224"
    python $PEM_PATH/../eiuss/eiuss.py gridsearch \
        --config            empirical_data/france_data/configs/config_${simulname}.py \
        --input             empirical_data/france_data/data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
        --outdir            empirical_data/france_data/${simulname} \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            1 \
        --operators         log10eta R0 \
        --range             " -8:-1.99:0.2;${2}:${3}:0.1" \
        --seed              202302${4}\
        ${5}
    ;;

  multi0101)
    echo $1 "| seed = " ${4}
    simulname="france_multiple_te=200101"
    python $PEM_PATH/../eiuss/eiuss.py gridsearch \
        --config            empirical_data/france_data/configs/config_${simulname}.py \
        --input             empirical_data/france_data/data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
        --outdir            empirical_data/france_data/${simulname} \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            1 \
        --operators         log10eta R0 \
        --range             " -8:-1.99:0.2;${2}:${3}:0.1" \
        --seed              202302${4}\
        ${5}
    ;;


  multi0108)
    echo $1 "| seed = " ${4}
    simulname="france_multiple_te=200108"
    python $PEM_PATH/../eiuss/eiuss.py gridsearch \
        --config            empirical_data/france_data/configs/config_${simulname}.py \
        --input             empirical_data/france_data/data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
        --outdir            empirical_data/france_data/${simulname} \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            1 \
        --operators         log10eta R0 \
        --range             " -8:-1.99:0.2;${2}:${3}:0.1" \
        --seed              202302${4}\
        ${5}

        ## range has blank at the start to provide - sign is working as a flag
    ;;

  esac

else      ## run from local
  case $1 in
  ready)
    sh ../run_remote.sh   bio_ready empirical_data/france_data/run_2D-loglikelihood
    sh ../run_remote.sh  rsph_ready empirical_data/france_data/run_2D-loglikelihood
    ;;

  run)
    #sh run_remote.sh   bio_run  session_2023-02-22.log   'for i in {0..9}   ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  single 1.0 2.5  ${i} --continued; done'
    #sh run_remote.sh   bio_run  session_2023-02-22.log   'for i in {10..19} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  single 2.5 4.0  ${i} --continued; done'
    #sh run_remote.sh   bio_run  session_2023-02-22.log   'for i in {20..29} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  single 4.0 4.51 ${i} --continued; done'

    #sh run_remote.sh   rsph_run session_2023-02-22.log   'for i in {0..9} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi1224 1.0 2.5 ${i} --continued; done'
    #sh run_remote.sh   rsph_run session_2023-02-22.log   'for i in {0..9} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi0101 1.0 2.5 ${i} --continued; done'
    #sh run_remote.sh   rsph_run session_2023-02-22.log   'for i in {0..9} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi0108 1.0 2.5 ${i} --continued; done'
    #sh run_remote.sh   rsph_run session_2023-02-22.log   'for i in {10..19} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi1224 2.5 4.01 ${i} --continued; done'
    #sh run_remote.sh   rsph_run session_2023-02-22.log   'for i in {10..19} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi0101 2.5 4.01 ${i} --continued; done'
    #sh run_remote.sh   rsph_run session_2023-02-22.log   'for i in {10..19} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi0108 2.5 4.01 ${i} --continued; done'

    sh ../run_remote.sh   rsph_run session_2023-02-22.log   'for i in {20..29} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi1224 4.1 5.0 ${i}; done'
    sh ../run_remote.sh   rsph_run session_2023-02-22.log   'for i in {20..29} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi0101 4.1 5.0 ${i}; done'
    sh ../run_remote.sh   rsph_run session_2023-02-22.log   'for i in {20..29} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi0108 4.1 5.0 ${i}; done'




    #sh run_remote.sh   bio_run session_2023-02-20.log   'for i in {0..9} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh   multi1224 0 ${i}; done'
    #sh run_remote.sh   bio_run session_2023-02-20.log   'for i in {0..9} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh   multi0101 0 ${i}; done'
    #sh run_remote.sh   rsph_run session_2023-02-20.log   'for i in {0..9} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi0108 0 ${i}; done'
    #sh run_remote.sh   rsph_run session_2023-02-20.log   'for i in {10..19} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi1224 1 ${i}; done'
    #sh run_remote.sh   rsph_run session_2023-02-20.log   'for i in {10..19} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi0101 1 ${i}; done'
    #sh run_remote.sh   rsph_run session_2023-02-20.log   'for i in {10..19} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi0108 1 ${i}; done'

    #sh run_remote.sh   rsph_run session_2023-02-20.log   'for i in {10 11 12 14 15 18 19} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi0108 1 ${i}; done'
    #sh run_remote.sh   rsph_run session_2023-02-20.log   'for i in 10 19;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi0108 1 ${i}; done'

    ;;

  download)
    echo "-------------";
    echo "rsph"
    sh ../run_remote.sh  rsph_download   empirical_data/france_data/france_multiple_te=200108/gridsearch_log10eta,R0
    sh ../run_remote.sh  rsph_download   empirical_data/france_data/france_multiple_te=200101/gridsearch_log10eta,R0
    sh ../run_remote.sh  rsph_download   empirical_data/france_data/france_multiple_te=191224/gridsearch_log10eta,R0

    echo "-------------";
    echo "bio"
    sh ../run_remote.sh   bio_download   empirical_data/france_data/france_single/gridsearch_timestart,R0
    ;;

  test)
    echo $1 "| seed = " ${4}
    simulname="france_multiple_te=200108"
    python $PEM_PATH/../../eiuss/eiuss.py gridsearch \
        --config            empirical_data/france_data/configs/config_${simulname}.py \
        --input             empirical_data/france_data/data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
        --outdir            empirical_data/france_data/${simulname} \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            1 \
        --operators         log10eta R0\
        --range             " -8:-1.99:0.2;${2}:${3}:0.1" \
        --seed              202302${4}\
        ${5}
      ;;

  esac

fi

