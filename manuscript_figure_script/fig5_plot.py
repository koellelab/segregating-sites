import os
import pickle
import datetime
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator

import density_plot
from fig2_plot import get_MLE_and_CI, n_segregating_over_time
from fig3_plot import read_and_get_mean_llk, plot_joint_estimation_95CI
# ---------------------------
def run(args):
    dir_simulated_dataset       = [args.simul_data_lowmu, args.simul_data_highmu]
    dir_llk_2Dgrid              = [args.llk_2dgrid_lowmu, args.llk_2dgrid_highmu]
    dir_phydyn_density_metadata = [args.phydyn_meta_lowmu, args.phydyn_meta_highmu]
    dir_phydyn_density_logFile  = [args.phydyn_log_lowmu, args.phydyn_log_highmu]

    # ---------------------------
    plt.rcParams.update({'font.size': 10})
    fig = plt.figure(figsize=[9, 8], constrained_layout=True, dpi = 300)
    gs0 = fig.add_gridspec(4, 2, width_ratios=[0.9, 1], height_ratios = [3,7, 3, 7], wspace=0.15, hspace = 0.01)

    gs00 = gs0[0, 0].subgridspec(1, 7, wspace=0.01)
    gs01 = gs0[1, 0].subgridspec(1, 1, wspace=0.01)
    gs02_1 = gs0[0, 1].subgridspec(7, 7, wspace=0.01)
    gs02_2 = gs0[1, 1].subgridspec(7, 7, wspace=0.01)

    b_gs00   = gs0[2, 0].subgridspec(1, 7, wspace=0.01)
    b_gs01   = gs0[3, 0].subgridspec(1, 1, wspace=0.01)
    b_gs02_1 = gs0[2, 1].subgridspec(7, 7, wspace=0.01)
    b_gs02_2 = gs0[3, 1].subgridspec(7, 7, wspace=0.01)

    axes1 = [fig.add_subplot(gs00[:, :6]),
             fig.add_subplot(gs01[:, :]),
             fig.add_subplot(gs02_1[2:, :5]),   ## upper
             fig.add_subplot(gs02_2[:, :5]),    ## center
             fig.add_subplot(gs02_2[:, 5:]),    ## right
             fig.add_subplot(gs02_1[:2, :5]),   ## upper
             fig.add_subplot(gs02_1[:2, 5:])    ## upper
            ]

    axes2 = [fig.add_subplot(b_gs00[:, :6]),
             fig.add_subplot(b_gs01[:, :]),
             fig.add_subplot(b_gs02_1[2:, :5]),   ## upper
             fig.add_subplot(b_gs02_2[:, :5]),    ## center
             fig.add_subplot(b_gs02_2[:, 5:]),    ## right
             fig.add_subplot(b_gs02_1[:2, :5]),   ## upper
             fig.add_subplot(b_gs02_1[:2, 5:])    ## upper
            ]


    # ---------------------------
    txt = open(f'manuscript figures/txt/{args.figname}.txt', 'w')

    for i, axes in enumerate([axes1, axes2]):
        # ---------------------------
        ## Panel A and D: segregating sites
        df_simulated_dataset = pd.read_csv(dir_simulated_dataset[i], sep="\t")
        n_segregating_over_time(axes[0], data=df_simulated_dataset)

        # ---------------------------
        ## Panel B and E
        df_2Dgrid_mean = read_and_get_mean_llk(dir_llk_2Dgrid[i], args.params, args.n_iter_per_cell)
        #df_2Dgrid_mean = df_2Dgrid_mean[df_2Dgrid_mean.timestart % 2 == 0]
        #df_2Dgrid_mean = df_2Dgrid_mean[df_2Dgrid_mean.timestart < 36]
        #df_2Dgrid_mean = df_2Dgrid_mean[df_2Dgrid_mean.R0 < 6.0 + 1e-8]

        df_2Dgrid_mean.mean_logL = np.where(~np.isinf(df_2Dgrid_mean.mean_logL), df_2Dgrid_mean.mean_logL, 1e5)
        df_2Dgrid_mean.mean_logL = np.where(~np.isnan(df_2Dgrid_mean.mean_logL), df_2Dgrid_mean.mean_logL, 1e5)
        df_2Dgrid_mean.mean_logL = np.where(df_2Dgrid_mean.sim_no == args.n_iter_per_cell, df_2Dgrid_mean.mean_logL, np.nan)





        ## heatmap for joint estimation
        plot_joint_estimation_95CI(axes[1], fig, df_2Dgrid_mean, args.params, args.true_params, gamma=25, majortick=10, df = 2)
        axes[1].set_xlabel(r'$R_0$', fontsize=14)
        axes[1].set_ylabel(r'$t_0$', fontsize=14)

        # print MLE and CI to txt
        df_2Dgrid_mean.mean_logL = np.where(df_2Dgrid_mean.mean_logL != 1e5, df_2Dgrid_mean.mean_logL, np.nan)

        mle_and_CI = get_MLE_and_CI(df_2Dgrid_mean, 'R0', MA=False, df=2)
        print("--- R0 ----\n", file=txt)
        print(mle_and_CI[0], "\n", file=txt)
        print(mle_and_CI[1], file=txt)
        print(mle_and_CI[2], "\n", file=txt)

        print(f"{mle_and_CI[0]['R0']:.3f}, {mle_and_CI[1]:.3f} - {mle_and_CI[2]:.3f}", file=txt)

        # print MLE and CI to txt; timestart
        mle_and_CI = get_MLE_and_CI(df_2Dgrid_mean, 'timestart', MA=False, df=2)
        print("--- timestart ----\n", file=txt)
        print(mle_and_CI[0], "\n", file=txt)
        print(mle_and_CI[1], file=txt)
        print(mle_and_CI[2], "\n", file=txt)

        print(f"{mle_and_CI[0]['timestart']:.3f}, {mle_and_CI[1]:.3f} - {mle_and_CI[2]:.3f}", file=txt)


        # ---------------------------
        ## Panel C and F
        density_plot.run(axes[2:5], dir_phydyn_density_metadata[i], dir_phydyn_density_logFile[i])

        ## axis range
        axes[3].set_xlim([1.0, 6.0])
        ## axis label
        axes[3].xaxis.set_major_locator(FixedLocator(np.arange(1, 6.1, 1)))  # FixedLocator
        axes[3].xaxis.set_minor_locator(FixedLocator(np.arange(1, 6.1, 0.2)))  # FixedLocator


        axes[3].set_ylim([-60, 38])
        axes[3].yaxis.set_major_locator(FixedLocator(np.arange(-60, 36, 10)))  # FixedLocator
        axes[3].yaxis.set_minor_locator(FixedLocator(np.arange(-60, 36, 2)))  # FixedLocator

        axes[2].set_xlim(axes[3].get_xlim())
        axes[4].set_ylim(axes[3].get_ylim())

    txt.close()
    # ---------------------------
    ## add subplot label
    labels = 'ABC DEF '.lower()
    for i, label in enumerate(labels):
        if i == 0:
            xloc, yloc, ax = (-0.35, 1.05, axes1[0])
        if i == 1:
            xloc, yloc, ax = (-0.2, 1.05, axes1[1])
        if i == 2:
            xloc, yloc, ax = (-0.35, 1.05, axes1[5])
        if i == 3:
            xloc, yloc, ax = (-0.35, 1.05, axes1[3])
        if i == 4:
            xloc, yloc, ax = (-0.35, 1.05, axes2[0])
        if i == 5:
            xloc, yloc, ax = (-0.2, 1.05, axes2[1])
        if i == 6:
            xloc, yloc, ax = (-0.35, 1.05, axes2[5])
        if i == 7:
            xloc, yloc, ax = (-0.35, 1.05, axes2[3])

        ax.text(xloc, yloc, label, transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')

    ## axis
    axes1[0].set_xlabel("Days post index case")
    axes1[0].set_ylabel("Number of\nsegregating\nsites")
    axes1[0].set_xticks(np.arange(52, 32, -4)[::-1])
    axes1[0].set_yticks([0, 5, 10])
    axes1[0].set_ylim([-0.3, 12])

    axes2[0].set_xlabel("Days post index case")
    axes2[0].set_ylabel("Number of\nsegregating\nsites")
    axes2[0].set_xticks(np.arange(52, 32, -4)[::-1])
    axes2[0].set_yticks([0, 5, 10, 15, 20])
    axes2[0].set_ylim([-0.3, 23])

    axes1[5].set_xticks([])
    axes1[5].set_yticks([])
    axes1[5].set_xlabel(None)
    axes1[5].set_ylabel(None)
    [axes1[5].spines[i].set_visible(False) for i in ['bottom', 'top', 'left', 'right']]

    axes1[6].set_xticks([])
    axes1[6].set_yticks([])
    axes1[6].set_xlabel(None)
    axes1[6].set_ylabel(None)
    [axes1[6].spines[i].set_visible(False) for i in ['bottom', 'top', 'left', 'right']]


    axes2[5].set_xticks([])
    axes2[5].set_yticks([])
    axes2[5].set_xlabel(None)
    axes2[5].set_ylabel(None)
    [axes2[5].spines[i].set_visible(False) for i in ['bottom', 'top', 'left', 'right']]

    axes2[6].set_xticks([])
    axes2[6].set_yticks([])
    axes2[6].set_xlabel(None)
    axes2[6].set_ylabel(None)
    [axes2[6].spines[i].set_visible(False) for i in ['bottom', 'top', 'left', 'right']]


    ## save figure
    fig.savefig(f'manuscript figures/pdf/{args.figname}.pdf')
    fig.savefig(f'manuscript figures/png/{args.figname}.png')

    fig.clf()
    plt.cla()

    print("done")
    return df_2Dgrid_mean








if __name__ == "__main__":
    import os
    import argparse


    if os.getcwd().find("manuscript_figure_script") == -1:
        parser = argparse.ArgumentParser()
        parser.add_argument('--figname', type=str, default="figure5")
        parser.add_argument('--params', type = str, nargs='+', default = ['R0', 'timestart'])
        parser.add_argument('--true_params', type = float, nargs='+')
        parser.add_argument('--n_iter_per_cell', type = int)

        parser.add_argument('--simul_data_lowmu', type=str, required=True)
        parser.add_argument('--simul_data_highmu', type=str, required=True)
        parser.add_argument('--llk_2dgrid_lowmu', type=str, required=True)
        parser.add_argument('--llk_2dgrid_highmu', type=str, required=True)

        parser.add_argument('--phydyn_meta_lowmu', type=str, required=False, default=None)
        parser.add_argument('--phydyn_meta_highmu', type=str, required=False, default=None)
        parser.add_argument('--phydyn_log_lowmu',   type=str, required=False, default=None)
        parser.add_argument('--phydyn_log_highmu',  type=str, required=False, default=None)


        args = parser.parse_args()

    print(args)
    hey = run(args)