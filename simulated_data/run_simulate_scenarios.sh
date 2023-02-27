#!/bin/bash
########################################
# 2023-01-27
# script for Figure 1
# simulate different scenarios and compare segregating sites trajectories
# : rng for mutations and epidemics are separated
#   - by setting param['rng2'] = param['rng'],
#     it works like an older version
########################################
PEM_PATH="$( cd "$( dirname "$0" )" && pwd -P )"
echo $PEM_PATH

# script below - - - - - - - - - - - - -
## Figure 1
case $1 in
  fig1)
    if false;then
    # --------------------------
    # figure 1AB: full, dense, sparse sampling
    simulname="simpleSEIR_R0=16e-1_mu2e-1"; seed_epi=1234
    mkdir -p simulated_data/${simulname}
    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       full  \
        --window_dt         4  \
        --seed              ${seed_epi}001 \
        --n_trial           1

    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       unif  \
        --sample_size       40  \
        --window_dt         4  \
        --seed              ${seed_epi}002 \
        --n_trial           30

    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       unif  \
        --sample_size       20  \
        --window_dt         4  \
        --seed              ${seed_epi}003 \
        --n_trial           30

    # --------------------------
    ### Fig 1C,D: higher R0
    simulname="simpleSEIR_R0=20e-1_mu2e-1"; seed_epi=1235
    mkdir -p simulated_data/${simulname}
    python $PEM_PATH/../eiuss/eiuss.py simulate_epidemic \
        --type              simple  \
        --config            simulated_data/configs/config_${simulname}.py \
        --seed              ${seed_epi}

    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       unif  \
        --sample_size       40  \
        --window_dt         4  \
        --seed              ${seed_epi}001 \
        --n_trial           30

    # --------------------------
    ### Fig 1E,F: transmission heterogeneity
    simulname="heteroSEIR_R0=16e-1_mu2e-1"; seed_epi=1236
    mkdir -p simulated_data/${simulname}
    python $PEM_PATH/../eiuss/eiuss.py simulate_epidemic \
        --type              simple  \
        --config            simulated_data/configs/config_${simulname}.py \
        --seed              ${seed_epi}

    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       unif  \
        --sample_size       40  \
        --window_dt         4  \
        --seed              ${seed_epi}001 \
        --n_trial           30

    fi
    # --------------------------
    ### Fig 1G,H: changing R0
    simulname="simpleSEIR_R0=16e-1_mu2e-1_reducedR0=75e-2"; seed_epi=1234
    mkdir -p simulated_data/${simulname}
    python $PEM_PATH/../eiuss/eiuss.py simulate_epidemic \
        --type              changingR0  \
        --config            simulated_data/configs/config_${simulname}.py \
        --seed              ${seed_epi}

    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       unif  \
        --sample_size       40  \
        --window_dt         4  \
        --seed              ${seed_epi}001 \
        --n_trial           30



    simulname="simpleSEIR_R0=16e-1_mu2e-1_reducedR0=11e-1"; seed_epi=1234
    mkdir -p simulated_data/${simulname}
    python $PEM_PATH/../eiuss/eiuss.py simulate_epidemic \
        --type              changingR0  \
        --config            simulated_data/configs/config_${simulname}.py \
        --seed              ${seed_epi}

    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       unif  \
        --sample_size       40  \
        --window_dt         4  \
        --seed              ${seed_epi}001 \
        --n_trial           30

    ;;

  figs1)
    simulname="simpleSEIR_R0=16e-1_mu2e-1"; seed_epi=1234
    mkdir -p simulated_data/${simulname}
    ## window 2, sparse
    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       unif  \
        --sample_size       10  \
        --window_dt         2  \
        --seed              ${seed_epi}004 \
        --n_trial           30

    ## window 2, dense
    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       unif  \
        --sample_size       20  \
        --window_dt         2  \
        --seed              ${seed_epi}005 \
        --n_trial           30

    ## window 6, sparse
    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       unif  \
        --sample_size       30  \
        --window_dt         6  \
        --seed              ${seed_epi}006 \
        --n_trial           30


    ## window 6, dense
    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       unif  \
        --sample_size       60  \
        --window_dt         6  \
        --seed              ${seed_epi}007 \
        --n_trial           30


    ;;

  figs2)
    simulname="heteroSEIR_R0=16e-1_mu2e-1_ph=3e-1"; seed_epi=1239
    mkdir -p simulated_data/${simulname}
    python $PEM_PATH/../eiuss/eiuss.py simulate_epidemic \
        --type              simple  \
        --config            simulated_data/configs/config_${simulname}.py \
        --seed              ${seed_epi}

    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       unif  \
        --sample_size       40  \
        --window_dt         4  \
        --seed              ${seed_epi}001 \
        --n_trial           30

    simulname="heteroSEIR_R0=16e-1_mu2e-1_ph=15e-2"; seed_epi=1239
    mkdir -p simulated_data/${simulname}
    python $PEM_PATH/../eiuss/eiuss.py simulate_epidemic \
        --type              simple  \
        --config            simulated_data/configs/config_${simulname}.py \
        --seed              ${seed_epi}

    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       unif  \
        --sample_size       40  \
        --window_dt         4  \
        --seed              ${seed_epi}001 \
        --n_trial           30

    ;;

  figs3)
  simulname="simpleSEIR_R0=16e-1_mu2e-1_reducedR0=75e-2_threshold=1000"; seed_epi=1234
  mkdir -p simulated_data/${simulname}
  python $PEM_PATH/../eiuss/eiuss.py simulate_epidemic \
      --type              changingR0  \
      --config            simulated_data/configs/config_${simulname}.py \
      --seed              ${seed_epi}

  python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
      --input             simulated_data/${simulname}/seed${seed_epi}  \
      --config            simulated_data/configs/config_${simulname}.py  \
      --sample_type       unif  \
      --sample_size       40  \
      --window_dt         4  \
      --seed              ${seed_epi}002 \
      --n_trial           30

  simulname="simpleSEIR_R0=16e-1_mu2e-1_reducedR0=11e-1_threshold=1000"; seed_epi=1234
  mkdir -p simulated_data/${simulname}
  python $PEM_PATH/../eiuss/eiuss.py simulate_epidemic \
      --type              changingR0  \
      --config            simulated_data/configs/config_${simulname}.py \
      --seed              ${seed_epi}

  python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
      --input             simulated_data/${simulname}/seed${seed_epi}  \
      --config            simulated_data/configs/config_${simulname}.py  \
      --sample_type       unif  \
      --sample_size       40  \
      --window_dt         4  \
      --seed              ${seed_epi}001 \
      --n_trial           30
    ;;


  figs4)
    simulname="simpleSEIR_R0=16e-1_mu2e-1"; seed_epi1=1234
    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi1}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       prop  \
        --sample_size       100  \
        --window_dt         4  \
        --seed              ${seed_epi1}001 \
        --n_trial           1
    ;;


  figs5)
    simulname="simpleSEIR_R0=16e-1_mu2e-1"; seed_epi1=1234
    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi1}  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       unif  \
        --sample_size       13  \
        --window_dt         4  \
        --seed              ${seed_epi1}001 \
        --n_trial           1

    ;;

  figs6)
    simulname="simpleSEIR_R0=16e-1_mu2e-1"; seed_epi1=1234
    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi1}  \
        --input_samp        simulated_data/${simulname}/seed${seed_epi1}_prop500_win4/seed230201  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       reset_window  \
        --sample_size       -1  \
        --window_dt         1  \
        --seed              ${seed_epi1}001 \
        --n_trial           1

    simulname="simpleSEIR_R0=16e-1_mu2e-1"; seed_epi1=1234
    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi1}  \
        --input_samp        simulated_data/${simulname}/seed${seed_epi1}_prop500_win4/seed230201  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       reset_window  \
        --sample_size       -1  \
        --window_dt         2  \
        --seed              ${seed_epi1}001 \
        --n_trial           1

    simulname="simpleSEIR_R0=16e-1_mu2e-1"; seed_epi1=1234
    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi1}  \
        --input_samp        simulated_data/${simulname}/seed${seed_epi1}_prop500_win4/seed230201  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       reset_window  \
        --sample_size       -1  \
        --window_dt         6  \
        --seed              ${seed_epi1}001 \
        --n_trial           1

    simulname="simpleSEIR_R0=16e-1_mu2e-1"; seed_epi1=1234
    python $PEM_PATH/../eiuss/eiuss.py simulate_sampling \
        --input             simulated_data/${simulname}/seed${seed_epi1}  \
        --input_samp        simulated_data/${simulname}/seed${seed_epi1}_prop500_win4/seed230201  \
        --config            simulated_data/configs/config_${simulname}.py  \
        --sample_type       reset_window  \
        --sample_size       -1  \
        --window_dt         10  \
        --seed              ${seed_epi1}001 \
        --n_trial           1

    ;;




  esac




