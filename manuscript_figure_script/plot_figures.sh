#!/bin/bash

PEM_PATH="$( cd "$( dirname "$0" )" && pwd -P )"
echo $PEM_PATH
# script below - - - - - -  - - - - - - -
case $1 in
## ------------------------------
## main text figures
## ------------------------------
  fig1)
    python $PEM_PATH/fig1_plot.py \
        --figname             figure1 \
        --baseline_statevar   simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_statevar.tsv \
        --baseline_full       simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_full_win4/seed1234001_segsites.tsv \
        --baseline_dense      simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif40_win4/seed1234002_trial##_segsites.tsv \
        --baseline_sparse     simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif20_win4/seed1234003_trial##_segsites.tsv \
          \
        --highR0              simulated_data/simpleSEIR_R0=20e-1_mu2e-1/seed1235_statevar.tsv \
        --highR0_dense        simulated_data/simpleSEIR_R0=20e-1_mu2e-1/seed1235_unif40_win4/seed1235001_trial##_segsites.tsv \
          \
        --transhet            simulated_data/heteroSEIR_R0=16e-1_mu2e-1/seed1236_statevar.tsv \
        --transhet_dense      simulated_data/heteroSEIR_R0=16e-1_mu2e-1/seed1236_unif40_win4/seed1236001_trial##_segsites.tsv \
          \
        --reduced0p75         simulated_data/simpleSEIR_R0=16e-1_mu2e-1_reducedR0=75e-2/seed1234_statevar.tsv \
        --reduced0p75_dense   simulated_data/simpleSEIR_R0=16e-1_mu2e-1_reducedR0=75e-2/seed1234_unif40_win4/seed1234001_trial##_segsites.tsv \
        --reduced1p1          simulated_data/simpleSEIR_R0=16e-1_mu2e-1_reducedR0=11e-1/seed1234_statevar.tsv \
        --reduced1p1_dense    simulated_data/simpleSEIR_R0=16e-1_mu2e-1_reducedR0=11e-1/seed1234_unif40_win4/seed1234001_trial##_segsites.tsv

    ;;

  fig2)         ## figure 2 - 1D log-likelihood
    python $PEM_PATH/fig2_plot.py \
        --simul_data          simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_segsites.tsv \
        --llk_rand            simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_gridsearch_R0 \
        --true_param          1.6 \
        --xlim                1 2.5 \
        --figname             figure2_R0=16e-1_seed1234

    ;;

  fig3)
    python $PEM_PATH/fig3_plot.py \
      --llk_2dgrid          simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_gridsearch_R0,timestart \
      --figname             figure3_R0=16e-1_seed1234_full
    ;;

  fig4)
    python $PEM_PATH/fig4_plot.py \
        --true_sim            simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_simtaj.npz \
        --true_segregating    simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_segsites.tsv \
        --reconstructed       simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_reconstruct \
        --figname             figure4_R0=16e-1_seed1234_reconstructed

    ;;

  fig5)
  python $PEM_PATH/fig5_plot.py \
        --figname             figure5_R0=16e-1_seed1234 \
        --params              R0 timestart \
        --true_params         1.6 0 \
        --n_iter_per_cell     20 \
        \
        --simul_data_lowmu    simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/seed230201_segsites.tsv \
        --llk_2dgrid_lowmu    simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/seed230201_gridsearch_R0,timestart \
        --simul_data_highmu   simulated_data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_unif10_win4_32-52/seed230201_segsites.tsv \
        --llk_2dgrid_highmu   simulated_data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_unif10_win4_32-52/seed230201_gridsearch_R0,timestart\
        \
        --phydyn_log_lowmu    simulated_data/phydyn_analysis/data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/simpleSEIR_R0=16e-1_mu2e-1_seed1234_unif10_win4_32-52.log \
        --phydyn_meta_lowmu   simulated_data/phydyn_analysis/data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/seed230201_sampling.tsv \
        --phydyn_log_highmu   simulated_data/phydyn_analysis/data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_unif10_win4_32-52/simpleSEIR_R0=16e-1_mu4e-1_seed221110_unif10_win4_32-52.log \
        --phydyn_meta_highmu  simulated_data/phydyn_analysis/data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_unif10_win4_32-52/seed230201_sampling.tsv

    ;;

  fig6)


    ;;

  fig7)
    python $PEM_PATH/fig7_plot.py \
        --figname             figure7 \
        --llk_2dgrid          empirical_data/france_data/france_single/gridsearch_timestart,R0 \
        --params              timestart R0 \
        --n_iter_per_cell     10
    ;;


  fig8)
    python $PEM_PATH/fig8_plot.py \
        --figname             figure8 \
        --params              log10eta R0 \
        --llk_2dgrid          empirical_data/france_data/france_multiple_te=191224/gridsearch_log10eta,R0 \
                              empirical_data/france_data/france_multiple_te=200101/gridsearch_log10eta,R0 \
                              empirical_data/france_data/france_multiple_te=200108/gridsearch_log10eta,R0 \
        --n_iter_per_cell     10
    ;;

  fig9)
    python $PEM_PATH/fig9_plot.py \
        --figname             figure9 \
        --observed_data       empirical_data/france_data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_match_largest_s.tsv \
        --dir_particle        empirical_data/france_data/France_multiple_exp_te=191224/reconstruct_2206051 \
                              empirical_data/france_data/France_multiple_exp_te=200101/reconstruct_2206051 \
                              empirical_data/france_data/France_multiple_exp_te=200108/reconstruct_2206051
    ;;

## ------------------------------
## supplementary figures
## ------------------------------
  figS1)
    ## Figure S1 - different time window lengths
    python $PEM_PATH/figS1_plot.py \
        --figname             figureS1 \
        --win2_sparse         simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win2/seed1234004_trial##_segsites.tsv \
        --win2_dense          simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif20_win2/seed1234005_trial##_segsites.tsv \
        --win4_sparse         simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif20_win4/seed1234003_trial##_segsites.tsv \
        --win4_dense          simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif40_win4/seed1234002_trial##_segsites.tsv \
        --win6_sparse         simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif30_win6/seed1234006_trial##_segsites.tsv \
        --win6_dense          simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif60_win6/seed1234007_trial##_segsites.tsv \
    ;;

  figS2)      ## Figure S2 - transmission heterogeneity
    python $PEM_PATH/figS2_plot.py \
        --figname             figureS2 \
        --statevar_080        simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_statevar.tsv \
        --segsite_080         simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif40_win4/seed1234002_trial##_segsites.tsv \
        --statevar_006        simulated_data/heteroSEIR_R0=16e-1_mu2e-1/seed1236_statevar.tsv \
        --segsite_006         simulated_data/heteroSEIR_R0=16e-1_mu2e-1/seed1236_unif40_win4/seed1236001_trial##_segsites.tsv \
        --statevar_015        simulated_data/heteroSEIR_R0=16e-1_mu2e-1_ph=15e-2/seed1239_statevar.tsv  \
        --segsite_015         simulated_data/heteroSEIR_R0=16e-1_mu2e-1_ph=15e-2/seed1239_unif40_win4/seed1239001_trial##_segsites.tsv \
        --statevar_030        simulated_data/heteroSEIR_R0=16e-1_mu2e-1_ph=3e-1/seed1239_statevar.tsv \
        --segsite_030         simulated_data/heteroSEIR_R0=16e-1_mu2e-1_ph=3e-1/seed1239_unif40_win4/seed1239001_trial##_segsites.tsv  \
    ;;

  figS3)      ## Figure S3 - threshold at 1000
    python $PEM_PATH/figS3_plot.py \
        --figname               figureS3  \
        --no_change_statevar    simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_statevar.tsv \
        --no_change_segsites    simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif40_win4/seed1234002_trial##_segsites.tsv \
        --change1p1_statevar    simulated_data/simpleSEIR_R0=16e-1_mu2e-1_reducedR0=11e-1_threshold=1000/seed1234_statevar.tsv \
        --change1p1_segsites    simulated_data/simpleSEIR_R0=16e-1_mu2e-1_reducedR0=11e-1_threshold=1000/seed1234_unif40_win4/seed1234001_trial##_segsites.tsv \
        --change0p75_statevar   simulated_data/simpleSEIR_R0=16e-1_mu2e-1_reducedR0=75e-2_threshold=1000/seed1234_statevar.tsv \
        --change0p75_segsites   simulated_data/simpleSEIR_R0=16e-1_mu2e-1_reducedR0=75e-2_threshold=1000/seed1234_unif40_win4/seed1234002_trial##_segsites.tsv \
    ;;

  figS4)    ## 1D log-likelhood with prop100
    python $PEM_PATH/figS4_plot.py \
        --simul_data          simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop100_win4/seed1234001_segsites.tsv \
        --llk_rand            simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop100_win4/seed1234001_gridsearch_R0 \
        --true_param          1.6 \
        --xlim                1 2.5 \
        --figname             figureS4_R0=16e-1_seed1234
    ;;

  figS5)  ## 1D log-likelhood with unif
    python $PEM_PATH/figS4_plot.py \
        --simul_data          simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif13_win4/seed1234001_segsites.tsv \
        --llk_rand            simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif13_win4/seed1234001_gridsearch_R0 \
        --true_param          1.6 \
        --xlim                1 2.5 \
        --figname             figureS5_R0=16e-1_seed1234

    ;;

  figS6)
    python $PEM_PATH/fig2_plot.py \
        --simul_data          simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif13_win4/seed1234001_segsites.tsv \
        --llk_rand            simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif13_win4/seed1234001_gridsearch_R0 \
        --true_param          1.6 \
        --xlim                1 2.5 \
        --figname             figureS6_R0=16e-1_seed1234


    ;;

  figS7)
    python $PEM_PATH/figS7_plot.py \
        --particle_data       simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_showparticles_R0,timestart/seed202302.tsv \
        --true_sim            simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_simtaj.npz \
        --true_data           simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_segsites.tsv \
        --figname             figureS7_R0=17e-1_timestart=16
    ;;


  figS8)
    python $PEM_PATH/figS8_plot.py \
        --figname       figureS8 \
        --dir_default   simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_simtaj.npz \
        --dir_high      simulated_data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_statevar.tsv


    ;;

  figS9)
    ## copied from empirical_data/france_data/scripts/run_all.sh
    cd $PEM_PATH/../empirical_data/france_data

    FranceSeqs=data/gisaid_hcov-19_2021_04_29_16.fasta
    python3 scripts/plot_empirical_data.py \
        --sDat ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
        --metadata ${FranceSeqs%.fasta}_date_country.csv \
        --outName "$PEM_PATH/../manuscript figures/figureS9"

    ;;


  figS12)
    python $PEM_PATH/figS12_plot.py \
        --figname       figureS12 \
        --dir_data      empirical_data/france_data/data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s_unambig.tsv \
                        empirical_data/france_data/data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
                        empirical_data/france_data/data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s_none.tsv


    ;;

esac
