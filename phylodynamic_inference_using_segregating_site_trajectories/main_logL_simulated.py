# ------------------
## main_logL.py
# ------------------
import sys
import os
import time
import shutil
import numpy as np
from importlib import import_module


orginal_dir = os.getcwd()
dir = "../"
sys.path.append(dir + "_script")

import tool
from SEIR_simple import particle
from SEIR_simple import likelihood as loglikelihood

# ---------------------------
# setup for sampling
data_setting = {'n_clade': 1,
                'n_class': 3,
                'window_dt': 4}
params_to_log = ["R0", "timeintro"]
# ---------------------------
input_file = sys.argv[1]        
config_file = sys.argv[2]
n_simul = int(sys.argv[3])
R0_start = float(sys.argv[4]); R0_end = float(sys.argv[5]); R0_interval = float(sys.argv[6]);
t0_start = float(sys.argv[7]); t0_end = float(sys.argv[8]); t0_interval = float(sys.argv[9]);


if len(sys.argv) > 10:
    tail = sys.argv[10]
else:
    tail = ""

# import and copy config file
input_dir = orginal_dir + "/input_data/" + input_file
output_dir = orginal_dir + "/output_result/" + config_file.split("/")[-1].replace(".py", "").replace("conf_", "") + tail + "/R0 = (%.2f %.2f %.2f), t0 = (%.2f %.2f %.2f)" %(R0_start, R0_end, R0_interval, t0_start, t0_end, t0_interval)

print ("input_dir :", input_dir)
print ("output_dir:",output_dir)

os.makedirs(output_dir, exist_ok=True)
shutil.copy(config_file, output_dir)
config = import_module(config_file.replace(".py", "").replace("/", "."))

R0_list = np.arange(R0_start, R0_end, R0_interval)
t0_list = np.arange(t0_start, t0_end, t0_interval)
print ("%d X %d = %d points" %(R0_list.size, t0_list.size, R0_list.size * t0_list.size))
# ---------------------------
landscape_dir = output_dir + "/logL_landscape_R0_" + str(n_simul) + "_" + str(R0_start) + "-" + str(R0_end) + "_" + str(t0_start) + "-" + str(t0_end)+ ".txt"
landscape_file = open(landscape_dir, 'w')
print('\t'.join([str(x) for x in params_to_log + ["logL"]]), file=landscape_file, flush=True)

# ---------------------------
imported_data = tool.import_data_n_segregating(input_dir)
params_PMCMC = config.set_param_pmcmc(imported_data)
# ----------------
i = 0
start = time.time()
for R0 in R0_list:
    for tstart in t0_list:

        print("<< R0: ", R0, ", timeintro: ", tstart, "(", i, "/", R0_list.size * t0_list.size, ") >>")
        params_PMCMC['timeintro'] = tstart
        params_PMCMC['R0'] = R0

        particles = [particle.Particle(imported_data['n_clades'], params_PMCMC['n_class']) for x in range(params_PMCMC['n_SMC_particles'])]  # allocate particles
        new_logL, particles, abandon_time_intro = loglikelihood.calculate(particles, params_PMCMC, imported_data, config)  # , params_to_log, log_file)

        print('\t'.join([str(x) for x in [R0, tstart, new_logL]]), file=landscape_file, flush=True)
        i += 1


end = time.time()
print("Total time : ", end - start)



