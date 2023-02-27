import os
import datetime
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator
from fig1_plot import plot_data
# ---------------------------
def run(win2_sparse, win2_dense, win4_sparse, win4_dense, win6_sparse, win6_dense, figname):

    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(1, 3, figsize=[8, 2.5], dpi=300)

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    # ---------------------------
    fnames = [[win2_sparse, win2_dense],
              [win4_sparse, win4_dense],
              [win6_sparse, win6_dense]]

    for i, window_size in enumerate([2, 4, 6]):
        for j, sampling in enumerate(["sparse", "dense"]):
            plot_data(axes[i], None, colors[j], fnames[i][j], linewidth=1, alpha=0.15)
            plot_data(axes[i], None, "k", fnames[i][j], linewidth=0.5, alpha=1, n_trial=1)     ## representative trajectory
            axes[i].plot([np.nan, np.nan], c=colors[j], label=sampling, linewidth=1, alpha=0.7)

        axes[i].legend()#title=f"samples per {window_size}-day window")
        axes[i].xaxis.set_major_locator(FixedLocator(np.arange(20, 168.1, 24)))  # FixedLocator
        axes[i].xaxis.set_minor_locator(FixedLocator(np.arange(20, 168.1, 8)))  # FixedLocator
        axes[i].set_xlim([20, 168])
        axes[i].grid(which = "both")

    # axis label

    axes[0].set_ylabel("Number of\nsegregating sites")
    axes[0].set_xlabel('Days post index case')
    axes[1].set_xlabel('Days post index case')
    axes[2].set_xlabel('Days post index case')

    axes[0].yaxis.set_major_locator(FixedLocator(np.arange(0, 300, 25)))  # FixedLocator
    axes[1].yaxis.set_major_locator(FixedLocator(np.arange(0, 300, 25)))  # FixedLocator
    axes[2].yaxis.set_major_locator(FixedLocator(np.arange(0, 300, 25)))  # FixedLocator

    ## add subplot label
    labels = 'ABC'.lower()
    for i, label in enumerate(labels):
        if i == 0:
            axes[i].text(-0.22, 1.1, label, transform=axes[i].transAxes, fontsize=16, fontweight='bold', va='top',
                         ha='left')
            axes[i].set_ylim([0, 300])
        else:
            axes[i].text(-0.18, 1.1, label, transform=axes[i].transAxes, fontsize=16, fontweight='bold', va='top',
                         ha='left')
            axes[i].set_ylim([0, 300])

    ## save figure
    fig.tight_layout(pad=0.5)
    fig.savefig(f'manuscript figures/pdf/{figname}.pdf')
    fig.savefig(f'manuscript figures/png/{figname}.png')

    fig.clf()
    plt.cla()

    return



if __name__ == "__main__":
    import os
    import argparse


    if os.getcwd().find("manuscript_figure_script") == -1:
        parser = argparse.ArgumentParser()
        parser.add_argument('--figname', type=str, default="figureS1")
        parser.add_argument('--win2_sparse', type=str)
        parser.add_argument('--win2_dense' , type=str)
        parser.add_argument('--win4_sparse', type=str)
        parser.add_argument('--win4_dense' , type=str)
        parser.add_argument('--win6_sparse', type=str)
        parser.add_argument('--win6_dense' , type=str)

        args = parser.parse_args()




    else:
        os.chdir(os.getcwd().replace("manuscript_figure_script", ""))
        parser = argparse.ArgumentParser()
        args = parser.parse_args()
        args.figname = "figureS1"

        args.win2_sparse = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/!!!!!!!_trial##_segsites.tsv"
        args.win2_dense  = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/!!!!!!!_trial##_segsites.tsv"
        args.win4_sparse = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/!!!!!!!_trial##_segsites.tsv"
        args.win4_dense  = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/!!!!!!!_trial##_segsites.tsv"
        args.win6_sparse = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/!!!!!!!_trial##_segsites.tsv"
        args.win6_dense  = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/!!!!!!!_trial##_segsites.tsv"





    print(args)
    run(args.win2_sparse,
        args.win2_dense,
        args.win4_sparse,
        args.win4_dense,
        args.win6_sparse,
        args.win6_dense,
        args.figname)
