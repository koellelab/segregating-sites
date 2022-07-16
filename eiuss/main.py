home_dir    = "/Users/yeongseon/Dropbox/PhD_Emory_university/2021-2022/Project_segregating_site"

import os
import sys
import pickle
import paramiko
import subprocess
import numpy as np
import pandas as pd
from importlib import import_module
# ---------------------------
nickname    = "test_github"
run_piuss          = True
plot_figure        = False
seed        = 220530
figure_date = 220530

working_dir = f'{nickname}'
os.makedirs(working_dir, exist_ok=True)
header      = f'{os.path.dirname(os.path.abspath(__file__))}/results/{__file__.split("/")[-1].replace(".py", "")}'
# ---------------------------
if run_piuss:

    f_dataset = f'{working_dir}/SEIR_simple_prop500_0928_seed123.tsv'
    f_config  = f'{working_dir}/config_test_github.py'
    f_out     = f'{working_dir}/{nickname}'
    f_llk     = f'{working_dir}/[mean_logL] {nickname}.tsv'

    ## copy config file to config folder
    import shutil
    shutil.copyfile(f'{home_dir}/{f_config}',
                    f'{home_dir}/piuss/config/{f_config.split("/")[-1]}')

    ## move to piuss script folder
    os.chdir(f"{home_dir}/piuss/");
    print(os.getcwd())
    my_env = os.environ.copy()

    ## run simulations through piuss
    for i in range(1):
        my_env["TAIL"]    = str(i)
        my_env["OUTFILE"] = f'{home_dir}/{f_out}'
        my_env["LLKFILE"] = f'{home_dir}/{f_llk}'

        run_simul = subprocess.Popen(
            [f"python", f"piuss.py",
             f"--inputData", f"../{f_dataset}",
             f"--config",    f"config/{f_config.split('/')[-1]}"
             ], env=my_env)  # , stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        output, errors = run_simul.communicate()
        run_simul.wait()

    os.chdir(f'{home_dir}/{working_dir}')
