#!/bin/bash
########################################
# 2023-01-27
# script for Figure 2
# calculate the likelihood for a range of R0
# : rng for mutations and epidemics are separated
#   - by setting param['rng2'] = param['rng'],
#     it works like an older version
########################################
PEM_PATH="$( cd "$( dirname "$0" )" && pwd -P )"
echo $PEM_PATH

# script below - - - - - -  - - - - - - -
if [ $(pwd | grep ypar279 | wc -l) -eq 1 ];then   ## run from server
  case $1 in
    run)
      simulname=$2        #"simpleSEIR_R0=16e-1_mu2e-1"
      dataname=$3         #"seed1234_prop500_win4/seed230201"
      ## full dataset with mu = 0.2
      python $PEM_PATH/eiuss/eiuss.py gridsearch \
          --config            simulated_data/configs/config_${simulname}.py \
          --input             simulated_data/${simulname}/${dataname}_segsites.tsv \
          --n_SMC_particles   200 \
          --n_grab            20 \
          --n_reps            1 \
          --operators         R0\
          --range             "1.0:2.51:0.05"\
          --seed              202302${5}
    ;;

    run2)   ## smaller bins size for range (1.4, 1.8)
      simulname=$2          #"simpleSEIR_R0=16e-1_mu2e-1"
      dataname=$3           #"seed1234_prop500_win4/seed230201"
      python $PEM_PATH/eiuss/eiuss.py gridsearch \
          --config            simulated_data/configs/config_${simulname}.py \
          --input             simulated_data/${simulname}/${dataname}_segsites.tsv \
          --n_SMC_particles   200 \
          --n_grab            20 \
          --n_reps            1 \
          --operators         R0\
          --range             "1.4:1.8:0.01"\
          --seed              202302${5}

  esac

else          ## run from local
  case $1 in
    ready)
      sh ../run_remote.sh   bio_ready simulated_data/run_1D-loglikelihood
      sh ../run_remote.sh  rsph_ready simulated_data/run_1D-loglikelihood
      ;;

    run)
      #sh run_remote.sh  bio_run        session_2023-02-12.log 'for i in {0..19}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh ${i} ${i}; done'
      #sh run_remote.sh  bio_run        session_2023-02-16.log 'for i in {20..24}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop500_win4/seed230201" 1 ${i}; done'
      #sh run_remote.sh  bio_run        session_2023-02-16.log 'for i in {25..29}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop500_win4/seed230201" 2 ${i}; done'
      #sh run_remote.sh  bio_run        session_2023-02-16.log 'for i in {30..34}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop500_win4/seed230201" 3 ${i}; done'
      #sh run_remote.sh  bio_run        session_2023-02-16.log 'for i in {35..39}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop500_win4/seed230201" 4 ${i}; done'

      #sh run_remote.sh  rsph_run        session_2023-02-12.log 'for i in {0..19}   ;do sbatch simulated_data/run_1D-loglikelihood_server.sh  run "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop100_win4/seed1234001" - ${i}; done'
      #sh run_remote.sh  rsph_run        session_2023-02-15.log 'for i in {20..24}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop100_win4/seed1234001" 1 ${i}; done'
      #sh run_remote.sh  rsph_run        session_2023-02-15.log 'for i in {25..29}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop100_win4/seed1234001" 2 ${i}; done'
      #sh run_remote.sh  rsph_run        session_2023-02-15.log 'for i in {30..34}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop100_win4/seed1234001" 3 ${i}; done'
      #sh run_remote.sh  rsph_run        session_2023-02-15.log 'for i in {35..39}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop100_win4/seed1234001" 4 ${i}; done'

      #sh run_remote.sh  rsph_run        session_2023-02-12.log 'for i in {0..19}   ;do sbatch simulated_data/run_1D-loglikelihood_server.sh  run "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_unif13_win4/seed1234001" - ${i}; done'
      #sh run_remote.sh  rsph_run        session_2023-02-15.log 'for i in {20..24}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_unif13_win4/seed1234001" 1 ${i}; done'
      #sh run_remote.sh  rsph_run        session_2023-02-15.log 'for i in {25..29}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_unif13_win4/seed1234001" 2 ${i}; done'
      #sh run_remote.sh  rsph_run        session_2023-02-15.log 'for i in {30..34}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_unif13_win4/seed1234001" 3 ${i}; done'
      #sh run_remote.sh  rsph_run        session_2023-02-15.log 'for i in {35..39}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_unif13_win4/seed1234001" 4 ${i}; done'


      ## for figure S6 ---------------------
      sh ../run_remote.sh  rsph_run ../session_2023-02-28.log 'for i in {0..19}   ;do sbatch simulated_data/run_1D-loglikelihood_server.sh   run "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop500_win4_resetwin1/seed230201/seed1234001" - ${i}; done'
      sh ../run_remote.sh  rsph_run ../session_2023-02-28.log 'for i in {0..19}   ;do sbatch simulated_data/run_1D-loglikelihood_server.sh   run "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop500_win4_resetwin2/seed230201/seed1234001" - ${i}; done'
      sh ../run_remote.sh  rsph_run ../session_2023-02-28.log 'for i in {0..19}   ;do sbatch simulated_data/run_1D-loglikelihood_server.sh   run "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop500_win4_resetwin6/seed230201/seed1234001" - ${i}; done'
      sh ../run_remote.sh  rsph_run ../session_2023-02-28.log 'for i in {0..19}   ;do sbatch simulated_data/run_1D-loglikelihood_server.sh   run "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop500_win4_resetwin10/seed230201/seed1234001" - ${i}; done'

      sh ../run_remote.sh  rsph_run ../session_2023-02-28.log 'for i in {20..39}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh   run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop500_win4_resetwin1/seed230201/seed1234001" - ${i}; done'
      sh ../run_remote.sh  rsph_run ../session_2023-02-28.log 'for i in {20..39}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh   run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop500_win4_resetwin2/seed230201/seed1234001" - ${i}; done'
      sh ../run_remote.sh  rsph_run ../session_2023-02-28.log 'for i in {20..39}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh   run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop500_win4_resetwin6/seed230201/seed1234001" - ${i}; done'
      sh ../run_remote.sh  rsph_run ../session_2023-02-28.log 'for i in {20..39}  ;do sbatch simulated_data/run_1D-loglikelihood_server.sh   run2 "simpleSEIR_R0=16e-1_mu2e-1" "seed1234_prop500_win4_resetwin10/seed230201/seed1234001" - ${i}; done'
      ;;


    download)
#      sh run_remote.sh   bio_download   simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_gridsearch_R0
#      sh run_remote.sh  rsph_download   simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop100_win4/seed1234001_gridsearch_R0
#      sh run_remote.sh  rsph_download   simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif13_win4/seed1234001_gridsearch_R0

      ## for figure S6 ---------------------
      sh ../run_remote.sh  rsph_download   simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4_resetwin1/seed230201
      sh ../run_remote.sh  rsph_download   simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4_resetwin2/seed230201
      sh ../run_remote.sh  rsph_download   simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4_resetwin6/seed230201
      sh ../run_remote.sh  rsph_download   simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4_resetwin10/seed230201
      ;;

    plot)
      python $PEM_PATH/../manuscript_figure_script/fig2_plot_grid.py \
          --simul_data          simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_segsites.tsv \
          --llk_rand            simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_gridsearch_R0 \
          --true_param          1.6 \
          --xlim                1 2.5 \
          --figname             figure2_R0=16e-1_seed1234
    ;;
  esac



fi

  ###### Done and uploaded to github repository #####