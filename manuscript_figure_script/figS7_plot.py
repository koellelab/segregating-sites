import os
import pickle
import datetime
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator

from fig4_plot import n_segregating_over_time, plot_state_variable
from fig7_plot import convert_to_matlab_date, convert_to_caldate
# ---------------------------









def run(args):
    txt = open(f'manuscript figures/txt/{args.figname}.txt', 'w')

    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(4, 1, figsize=[8, 8], dpi=300, sharex = True)
    # ---------------------------
    df_simulated_dataset = pd.read_csv(args.true_data, sep = "\t")
    df_true_statevar     = np.load    (args.true_sim, allow_pickle = True)["statevar"]
    df_particle          = pd.read_csv(args.particle_data, sep = "\t")

    color = {'X': 'grey', 'V': 'k'}
    alpha = {'X': 0.1, 'V': 0.3}
    ## (A, B, C, D) reconstructed state variables
    for i in range(len(df_particle)):
        for j, state in enumerate(['seg', 's', 'e', 'i']):
            axes[j].plot([df_particle.iloc[i]['t_before'], df_particle.iloc[i]['t_after']],
                      [df_particle.iloc[i][f'{state}_before'], df_particle.iloc[i][f'{state}_after']],
                      c=color[df_particle.iloc[i]['selected']],
                      alpha=alpha[df_particle.iloc[i]['selected']],
                      linewidth=1)



    ## true state variables
    n_segregating_over_time(axes[0], df_simulated_dataset, particle=False)
    plot_state_variable    (axes[1], df_true_statevar, 1, 'susceptible', particle=False)    #  [0, 110000])
    plot_state_variable    (axes[2], df_true_statevar, 2, 'exposed',  particle=False)        #  [0, 6000])
    plot_state_variable    (axes[3], df_true_statevar, 3, 'infected',  particle=False)       #  [0, 8000])

    # ---------------------------
    labels = 'ABCD'.lower()
    for i, label_row in enumerate(labels):
        axes[i].text(-0.15, 1.1, label_row, transform=axes[i].transAxes, fontsize=16, fontweight='bold', va='top',
                     ha='left')

    ## save figure
    fig.tight_layout()
    fig.savefig(f'manuscript figures/pdf/{args.figname}.pdf')
    fig.savefig(f'manuscript figures/png/{args.figname}.png')

    fig.clf()
    plt.cla()

    return

if __name__ == "__main__":
    import os
    import argparse


    if os.getcwd().find("manuscript_figure_script") == -1:
        parser = argparse.ArgumentParser()
        parser.add_argument('--figname', type=str, default="figureS7")
        parser.add_argument('--true_sim', type=str)
        parser.add_argument('--true_data', type=str)
        parser.add_argument('--particle_data', type=str)


        args = parser.parse_args()

    else:
        os.chdir(os.getcwd().replace("manuscript_figure_script", ""))
        parser = argparse.ArgumentParser()
        args = parser.parse_args()

    print(args)
    run(args)
