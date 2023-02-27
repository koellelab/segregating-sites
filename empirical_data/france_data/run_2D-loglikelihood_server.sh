#!/bin/bash
#SBATCH --partition=3day-long
#SBATCH --output=slurm-%A.%a.out
#SBATCH --error=slurm-%A.%a.err
#SBATCH --job-name=randomsearch
#----------------
source /opt/anaconda/etc/profile.d/conda.sh
conda activate /apps/conda/ypar279/envs/project_ss

if [ -n $SLURM_JOB_ID ] ; then
THEPATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
else
THEPATH=$(realpath $0)
fi

PEM_PATH=$(dirname $(dirname "${THEPATH}"))
echo $PEM_PATH






# script below - - - - - -  - - - - - - -
if [ $(pwd | grep ypar279 | wc -l) -eq 1 ];then   ## run from server
  case $1 in
  single)
    echo $1
    ## full dataset with mu = 0.2
    simulname="france_single"
    python $PEM_PATH/../eiuss/eiuss.py gridsearch \
        --config            empirical_data/france_data/configs/config_${simulname}.py \
        --input             empirical_data/france_data/data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
        --outdir            empirical_data/france_data/${simulname} \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            1 \
        --operators         R0 timestart\
        --range             "1.0:4.51:0.1;2020/02/17:2019/12/01:-3" \
        --seed              202302${3}\
        ${4}
      ;;

  multi1224)
    echo $1
    simulname="france_multiple_te=191224"
    python $PEM_PATH/../eiuss/eiuss.py gridsearch \
        --config            empirical_data/france_data/configs/config_${simulname}.py \
        --input             empirical_data/france_data/data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
        --outdir            empirical_data/france_data/${simulname} \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            1 \
        --operators         R0 log10eta\
        --range             "1.${2}:4.01:0.2;-8:-1.99:0.2" \
        --seed              202302${3}\
        ${4}
    ;;

  multi0101)
    echo $1
    simulname="france_multiple_te=200101"
    python $PEM_PATH/../eiuss/eiuss.py gridsearch \
        --config            empirical_data/france_data/configs/config_${simulname}.py \
        --input             empirical_data/france_data/data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
        --outdir            empirical_data/france_data/${simulname} \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            1 \
        --operators         R0 log10eta\
        --range             "1.${2}:4.01:0.2;-8:-1.99:0.2" \
        --seed              202302${3}\
        ${4}
    ;;


  multi0108)
    echo $1
    simulname="france_multiple_te=200108"
    python $PEM_PATH/../eiuss/eiuss.py gridsearch \
        --config            empirical_data/france_data/configs/config_${simulname}.py \
        --input             empirical_data/france_data/data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
        --outdir            empirical_data/france_data/${simulname} \
        --n_SMC_particles   200 \
        --n_grab            20 \
        --n_reps            1 \
        --operators         R0 log10eta\
        --range             "1.${2}:4.01:0.2;-8:-1.99:0.2" \
        --seed              202302${3}\
        ${4}
    ;;

  esac

else      ## run from local
  case $1 in
  ready)
    sh run_remote.sh   bio_ready empirical_data/france_data/run_2D-loglikelihood
    sh run_remote.sh  rsph_ready empirical_data/france_data/run_2D-loglikelihood
    ;;

  run)
    #sh run_remote.sh   bio_run session_2023-02-20.log   'for i in {0..9} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  single _ ${i}; done'
    sh run_remote.sh   bio_run session_2023-02-20.log   'for i in {0..9} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh   multi1224 0 ${i}; done'
    sh run_remote.sh   bio_run session_2023-02-20.log   'for i in {0..9} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh   multi0101 0 ${i}; done'
    sh run_remote.sh   rsph_run session_2023-02-20.log   'for i in {0..9} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi0108 0 ${i}; done'
    sh run_remote.sh   rsph_run session_2023-02-20.log   'for i in {10..19} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi1224 1 ${i}; done'
    sh run_remote.sh   rsph_run session_2023-02-20.log   'for i in {10..19} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi0101 1 ${i}; done'
    sh run_remote.sh   rsph_run session_2023-02-20.log   'for i in {10..19} ;do sbatch empirical_data/france_data/run_2D-loglikelihood_server.sh  multi0108 1 ${i}; done'
    ;;

  download)
    sh run_remote.sh  rsph_download   empirical_data/france_data/
    sh run_remote.sh   bio_download   empirical_data/france_data/
    ;;

  test)
    echo $1
    ## full dataset with mu = 0.2
    simulname="france_multiple_te=200108"
    python $PEM_PATH/../../eiuss/eiuss.py gridsearch \
        --config            empirical_data/france_data/configs/config_${simulname}.py \
        --input             empirical_data/france_data/data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
        --outdir            empirical_data/france_data/${simulname} \
        --n_SMC_particles   1 \
        --n_grab            20 \
        --n_reps            1 \
        --operators         R0 log10eta\
        --range             "1.${2}:4.01:0.2;-8:-1.99:0.2" \
        --seed              202302${3}\
        ${4}
      ;;

  esac

fi


