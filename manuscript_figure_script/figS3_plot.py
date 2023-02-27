import os
import datetime
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator

from fig1_plot import plot_data, plot_statevar
# ---------------------------
def run(args):
    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(2, 1, figsize=[4, 5.5], dpi=300)

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    # ---------------------------
    plot_statevar(axes[0], r'R$_{0}$ unchanged',       "k",       args.no_change_statevar)
    plot_data    (axes[1], r'R$_{0}$ unchanged',       "k",       args.no_change_segsites)
    plot_statevar(axes[0], r'R$_{0}$ reduced to 1.1' , "#ff5f43", args.change1p1_statevar)
    plot_data    (axes[1], r'R$_{0}$ reduced to 1.1' , "#ff5f43", args.change1p1_segsites)
    plot_statevar(axes[0], r'R$_{0}$ reduced to 0.75', "#fcb708", args.change0p75_statevar)
    plot_data    (axes[1], r'R$_{0}$ reduced to 0.75', "#fcb708", args.change0p75_segsites)


    # ---------------------------
    axes[0].legend(loc="upper right")

    # xticks
    for i in [0, 1]:
        axes[i].xaxis.set_major_locator(FixedLocator(np.arange(00, 200.1, 24)))  # FixedLocator
        axes[i].xaxis.set_minor_locator(FixedLocator(np.arange(00, 200.1, 8)))  # FixedLocator
        axes[i].set_xlim([0, 200])
        axes[i].grid(which="both")

    # yticks
    axes[0].set_ylim([0, 10000]);
    axes[0].ticklabel_format(axis='y', style="sci", scilimits=(3, 4))
    axes[1].set_ylim([0, 225]);

    # axis label
    axes[0].set_ylabel("Number of infected individuals")
    axes[1].set_ylabel("Number of segregating sites")
    axes[1].set_xlabel('Days post index case')

    ## add subplot label
    labels = 'AB'.lower()
    for i, label in enumerate(labels):
        axes[i].text(-0.22, 1.1, label, transform=axes[i].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')




    ## save figure
    fig.tight_layout(pad=0.5)
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
        parser.add_argument('--figname', type=str, default="figureS3")
        parser.add_argument('--no_change_statevar', type=str)
        parser.add_argument('--no_change_segsites' , type=str)
        parser.add_argument('--change1p1_statevar', type=str)
        parser.add_argument('--change1p1_segsites' , type=str)
        parser.add_argument('--change0p75_statevar', type=str)
        parser.add_argument('--change0p75_segsites' , type=str)

        args = parser.parse_args()



    else:
        os.chdir(os.getcwd().replace("manuscript_figure_script", ""))
        parser = argparse.ArgumentParser()
        args = parser.parse_args()
        args.figname = "figureS3"
        args.no_change_statevar  = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed20230201_statevar.tsv"
        args.no_change_segsites  = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed20230201_unif40_win4/seed230201_trial##_segsites.tsv"
        args.change1p1_statevar  = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1_reducedR0=11e-1_threshold=1000/seed20230201_statevar.tsv"
        args.change1p1_segsites  = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1_reducedR0=11e-1_threshold=1000/seed20230201_unif40_win4/seed230201_trial##_segsites.tsv"
        args.change0p75_statevar = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1_reducedR0=75e-2_threshold=1000/seed20230201_statevar.tsv"
        args.change0p75_segsites = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1_reducedR0=75e-2_threshold=1000/seed20230201_unif40_win4/seed230201_trial##_segsites.tsv"



    print(args)
    run(args)