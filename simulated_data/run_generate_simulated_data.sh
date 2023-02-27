#!/bin/bash
########################################
# simulate epidemics and sampling
# : rng for mutations and epidemics are separated
#   - by setting param['rng2'] = param['rng'],
#     it works like an older version
########################################
PEM_PATH="$( cd "$( dirname "$0" )" && pwd -P )"
echo $PEM_PATH

# script below - - - - - -  - - - - - - -
## 1. simulate epidemiological dynamics
##    : this should reproduce the simulations in the manuscript
if false;then
simulname="simpleSEIR_R0=16e-1_mu2e-1"
  mkdir -p simulated_data/${simulname}
  ## simulate an epidemic
  seed_epi1=1234
  python $PEM_PATH/../eiuss/eiuss.py simulate_epidemic \
      --type              simple  \
      --config            simulated_data/configs/config_${simulname}.py \
      --seed              ${seed_epi1} #> simulated_data/${simulname}/simulate_epidemic_seed${seed_epi}.txt

  python $PEM_PATH/../eiuss/eiuss_test.py show_epidynamics \
    simulated_data/${simulname}/seed${seed_epi1}_statevar.tsv

simulname="simpleSEIR_R0=16e-1_mu4e-1"
  mkdir -p simulated_data/${simulname}
  ## simulate an epidemic
  seed_epi2=221110
  python $PEM_PATH/../eiuss/eiuss.py simulate_epidemic \
      --type              simple  \
      --config            simulated_data/configs/config_${simulname}.py \
      --seed              ${seed_epi2} #> simulated_data/${simulname}/simulate_epidemic_seed${seed_epi}.txt

  python $PEM_PATH/../eiuss/eiuss_test.py show_epidynamics \
    simulated_data/${simulname}/seed${seed_epi2}_statevar.tsv
fi



simulname="simpleSEIR_R0=16e-1_mu2e-1"
  mkdir -p simulated_data/${simulname}
  ## simulate an epidemic
  seed_epi1=1234
  python $PEM_PATH/../eiuss/eiuss.py simulate_epidemic \
      --type              simple  \
      --config            simulated_data/configs/config_${simulname}.py \
      --seed              ${seed_epi1} #> simulated_data/${simulname}/simulate_epidemic_seed${seed_epi}.txt

  python $PEM_PATH/../eiuss/eiuss_test.py show_epidynamics \
    simulated_data/${simulname}/seed${seed_epi1}_statevar.tsv







## 2. generate segregating site data
##    : different from the 2022-Dec version as we fixed the bug
##    : sampling was done with replacement while picking an individual
if false;then
simulname="simpleSEIR_R0=16e-1_mu2e-1"; seed_epi1=1234
  ## simulate sampling - "full dataset"
  seed_sam=230201
  python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
      --input             simulated_data/${simulname}/seed${seed_epi1}  \
      --config            simulated_data/configs/config_${simulname}.py  \
      --sample_type       prop  \
      --sample_size       500  \
      --window_dt         4  \
      --seed              ${seed_sam} \
      --n_trial           1   #> simulated_data/${simulname}/simulate_sampling_seed${seed_epi}.txt

  python $PEM_PATH/../eiuss/eiuss_test.py check_segsite_trajctories \
    simulated_data/${simulname}/seed${seed_epi1}_prop500_win4/seed${seed_sam}_segsites.tsv

  ## simulate sampling - "short dataset"
  seed_sam=230201
  python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
      --input             simulated_data/${simulname}/seed${seed_epi1}  \
      --config            simulated_data/configs/config_${simulname}.py  \
      --sample_type       unif  \
      --sample_between    32 52 \
      --sample_size       10  \
      --window_dt         4  \
      --seed              ${seed_sam} \
      --n_trial           1   #> simulated_data/${simulname}/simulate_sampling_seed${seed_epi}.txt

  python $PEM_PATH/../eiuss/eiuss_test.py check_segsite_trajctories \
    simulated_data/${simulname}/seed${seed_epi1}_unif10_win4_32-52/seed${seed_sam}_segsites.tsv


simulname="simpleSEIR_R0=16e-1_mu4e-1"; seed_epi2=221110
  seed_sam=230201
  python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
      --input             simulated_data/${simulname}/seed${seed_epi2}  \
      --config            simulated_data/configs/config_${simulname}.py  \
      --sample_type       unif  \
      --sample_between    32 52 \
      --sample_size       10  \
      --window_dt         4  \
      --seed              ${seed_sam} \
      --n_trial           1   #> simulated_data/${simulname}/simulate_sampling_seed${seed_epi}.txt

  python $PEM_PATH/../eiuss/eiuss_test.py check_segsite_trajctories \
    simulated_data/${simulname}/seed${seed_epi2}_unif10_win4_32-52/seed${seed_sam}_segsites.tsv

fi








