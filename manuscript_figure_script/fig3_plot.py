import os
import pickle
import datetime
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator

import density_plot
from fig2_plot import get_MLE_and_CI, read_and_get_mean_llk
# ---------------------------


def plot_joint_estimation_95CI(ax, fig,  df_surface, params, trueparam, gamma, tick=[], majortick = 3, majortick_y = 5, edgecolor = "w", alpha = 1, ax2 = None, df = 1, xlim =None, ylim = None):
    import matplotlib.colors as colors
    df_surface_ = df_surface
    df_surface  = df_surface.pivot(columns=params[0], index=params[1], values='mean_logL')

    # https://matplotlib.org/stable/tutorials/colors/colormapnorms.html
    # https: // rfriend.tistory.com / 419
    cmap = plt.get_cmap('viridis').copy()
    cmap.set_under(color='w', alpha=1.)
    cmap.set_over(color='w', alpha=1.)
    cmap.set_bad(color='lightgray', alpha=1.)

    '''
    p = ax.pcolor(df_surface, cmap = cmap,
                    norm=colors.PowerNorm(gamma=gamma,
                                          vmin=np.nanmin(df_surface.values[~np.isinf(df_surface.values)]),
                                          vmax=np.nanmax(df_surface.values)),
                    shading='auto', alpha = alpha,
                    edgecolor=edgecolor, linewidth=0.00005)  # https: // rfriend.tistory.com / 419
    '''


    if df == 1:
        thres = 1.92
    elif df == 2:
        thres = 2.9957325

    p = ax.pcolormesh(df_surface,
                      cmap=cmap,
                      norm=colors.PowerNorm(gamma=gamma,
                                            vmin=np.nanmin(df_surface.values[~np.isinf(df_surface.values)]),
                                            vmax=np.nanmax(df_surface.values[df_surface.values != 1e5])),
                       shading='auto', alpha = 1,linewidth=0.0025)
                       #edgecolor=edgecolor,


    ## color bars for the heatmap
    ticks = [np.nanmin(df_surface_[~np.isinf(df_surface_)].mean_logL),
             np.nanmax(df_surface_[df_surface_ != 1e5].mean_logL) - thres,
             np.nanmax(df_surface_[df_surface_ != 1e5].mean_logL)] + tick
    cbar = fig.colorbar(p, ax=ax, ticks=ticks)
    ax.set_aspect('auto')
    cbar.ax.tick_params(labelsize=8)


    ## highlight  the 95% confidence interval
    max_mean_logL = df_surface_[df_surface_.mean_logL != 1e5].mean_logL.max()

    from matplotlib.collections import LineCollection

    x = np.arange(0, len(df_surface.columns)+1, 1)
    y = np.arange(0, len(df_surface.index)+1, 1)
    z = df_surface.copy()
    z = np.where(z != 1e5, z, np.nan)

    def add_iso_line(ax, value, color):
        ## https://stackoverflow.com/questions/63458863/way-to-contour-outer-edge-of-selected-grid-region-in-python
        v = np.diff(z > value, axis=1)
        h = np.diff(z > value, axis=0)

        l = np.argwhere(v.T)
        vlines = np.array(list(zip(np.stack((x[l[:, 0] + 1], y[l[:, 1]])).T,
                                   np.stack((x[l[:, 0] + 1], y[l[:, 1] + 1])).T)))
        l = np.argwhere(h.T)
        l_ = l[~(l[:, 0] + 1 < len(x)) * (l[:, 1] + 1 < len(y))]
        l = l[(l[:, 0]+1 < len(x)) *  (l[:, 1]+1 < len(y))]
        hlines = np.array(list(zip(np.stack((x[l[:, 0]],     y[l[:, 1] + 1])).T,
                                   np.stack((x[l[:, 0] + 1], y[l[:, 1] + 1])).T)) +
                          list(zip(np.stack((x[l_[:, 0]], y[l_[:, 1] + 1])).T,
                                   np.stack((x[l_[:, 0]], y[l_[:, 1] + 1])).T)),)

        lines = np.vstack((vlines, hlines))
        ax.add_collection(LineCollection(lines, lw=1, colors=color))



    add_iso_line(ax, max_mean_logL-thres, 'r')

    """
    max_mean_logL   = df_surface_.mean_logL.max()
    df_surface_95CI = df_surface_.copy()
    df_surface_95CI = df_surface_95CI.pivot(columns=params[0], index=params[1], values='mean_logL')
    df_surface_95CI = np.where(df_surface_95CI > max_mean_logL - 2,
                               df_surface_95CI, np.nan)

    p = ax.pcolor(df_surface_95CI, alpha=1, edgecolor='r', linewidth=0.5)  # https: // rfriend.tistory.com / 419
    """

    ## first paramter on x axis
    ax.set_xlabel(params[0], fontsize=14)
    print(df_surface.columns)
    ax.xaxis.set_major_locator(FixedLocator(np.arange(0.5, len(df_surface.columns), 1)))  # FixedLocator
    ax.xaxis.set_minor_locator(FixedLocator(np.arange(0.5, len(df_surface.columns), 1)))  # FixedLocator
    ax.set_xticklabels(["%.1f" % (df_surface.columns[i]) if i % majortick == 0 else "" for i in range(len(df_surface.columns))])
    ax.xaxis.set_major_locator(FixedLocator(np.arange(0.5, len(df_surface.columns), majortick)))  # FixedLocator

    #ax.set_xticklabels(["%.1f" % (df_surface.columns[i]) for i in range(len(df_surface.columns)) if i % majortick == 0])
    #
    #ax.xaxis.set_major_locator(FixedLocator(np.arange(0.5, len(df_surface.columns), majortick)))  # FixedLocator

    ## true line
    if not np.isnan(trueparam[0]):
        loc = np.arange(0.5, len(df_surface.columns), 1)[np.where(df_surface.columns == trueparam[0])]
        ax.axvline(x=loc, linewidth=1.5, ls='--')

    ##set x and y range
    if not xlim is None:
        loc1 = np.arange(0.5, len(df_surface.columns), 1)[np.where(df_surface.columns == xlim[0])]
        loc2 = np.arange(0.5, len(df_surface.columns), 1)[np.where(df_surface.columns == xlim[1])]
        ax.set_xlim([loc1-0.5, loc2+0.5])
    #ax.set_xlim([0, len(df_surface.columns)])





    ## second paramter on y axis
    ax.set_ylabel(params[1], fontsize=14)
    print(df_surface.index)
    ax.yaxis.set_major_locator(FixedLocator(np.arange(0.5, len(df_surface.index), 1)))  # FixedLocator
    ax.yaxis.set_minor_locator(FixedLocator(np.arange(0.5, len(df_surface.index), 1)))  # FixedLocator
    ax.set_yticklabels(["%d" % (df_surface.index[i]) if i % majortick_y == 0 else "" for i in range(len(df_surface.index))])
    ax.yaxis.set_major_locator(FixedLocator(np.arange(0.5, len(df_surface.index), majortick_y)))  # FixedLocator



    ## true line
    if not np.isnan(trueparam[1]):
        loc = np.arange(0.5, len(df_surface.index), 1)[np.where(df_surface.index == trueparam[1])]
        ax.axhline(y=loc, linewidth=1.5, ls='--')

    ##set x and y range
    if not ylim is None:
        loc1 = np.arange(0.5, len(df_surface.index), 1)[np.where(df_surface.index == ylim[0])]
        loc2 = np.arange(0.5, len(df_surface.index), 1)[np.where(df_surface.index == ylim[1])]
        ax.set_ylim([loc1-0.5, loc2+0.5])



    return



def run(dir_llk_2Dgrid, dir_phydyn_density_metadata, dir_phydyn_density_logFile, figname):

    # ---------------------------
    plt.rcParams.update({'font.size': 10})
    fig = plt.figure(figsize=[9, 4], constrained_layout=True, dpi = 300)

    gs0 = fig.add_gridspec(1, 2, width_ratios=[0.9, 1], wspace=8)
    gs00 = gs0[0].subgridspec(7, 7, wspace=0.01)
    gs01 = gs0[1].subgridspec(7, 7, wspace=0.01)

    axes = [fig.add_subplot(gs00[2:, :]),
            fig.add_subplot(gs01[:2, :5]),
            fig.add_subplot(gs01[2:, :5]),
            fig.add_subplot(gs01[2:, 5:])]

    # ---------------------------
    ## panel A
    n_iter_per_cell = 20; df = 2
    df_2Dgrid_mean = read_and_get_mean_llk(dir_llk_2Dgrid, ['R0', 'timestart'], n_iter_per_cell)
    df_2Dgrid_mean.to_csv(dir_llk_2Dgrid + "/meanllk.tsv", sep = "\t", index = False)

    df_2Dgrid_mean.mean_logL = np.where(~np.isinf(df_2Dgrid_mean.mean_logL), df_2Dgrid_mean.mean_logL, 1e5)
    df_2Dgrid_mean.mean_logL = np.where(~np.isnan(df_2Dgrid_mean.mean_logL), df_2Dgrid_mean.mean_logL, 1e5)
    df_2Dgrid_mean.mean_logL = np.where(df_2Dgrid_mean.sim_no == n_iter_per_cell, df_2Dgrid_mean.mean_logL, np.nan)

    plot_joint_estimation_95CI(axes[0], fig, df_2Dgrid_mean, ['R0', 'timestart'], [1.6, 0], majortick = 5,gamma=12, tick=[-80, -90], df = df)
    
    ## axis
    axes[0].set_xlabel(r'$R_0$', fontsize=14)
    axes[0].set_ylabel(r'$t_0$', fontsize=14)

    # print MLE and CI to txt; R0
    df_2Dgrid_mean.mean_logL = np.where(df_2Dgrid_mean.mean_logL != 1e5, df_2Dgrid_mean.mean_logL, np.nan)
    
    txt = open(f'manuscript figures/txt/{figname}.txt', 'w')
    mle_and_CI = get_MLE_and_CI(df_2Dgrid_mean, 'R0', MA=False, df=2)
    print("--- R0 ----\n", file=txt)
    print(mle_and_CI[0], "\n",  file = txt)
    print(mle_and_CI[1], file = txt)
    print(mle_and_CI[2], "\n", file = txt)

    print(f"{mle_and_CI[0]['R0']:.3f}, {mle_and_CI[1]:.3f} - {mle_and_CI[2]:.3f}", file = txt)

    # print MLE and CI to txt; timestart
    mle_and_CI = get_MLE_and_CI(df_2Dgrid_mean, 'timestart', MA=False, df = 2)
    print("--- timestart ----\n", file=txt)
    print(mle_and_CI[0], "\n", file=txt)
    print(mle_and_CI[1], file=txt)
    print(mle_and_CI[2], "\n", file=txt)

    print(f"{mle_and_CI[0]['timestart']:.3f}, {mle_and_CI[1]:.3f} - {mle_and_CI[2]:.3f}", file=txt)

    txt.close()


    # --------------------------------
    ## panel B
    if dir_phydyn_density_metadata and dir_phydyn_density_logFile:
        density_plot.run(axes[1:4], dir_phydyn_density_metadata, dir_phydyn_density_logFile)

    ## axis range
    axes[2].set_xlim([1.0, 2.5])
    axes[2].set_ylim([-50, 38])
    axes[1].set_xlim(axes[2].get_xlim())
    axes[3].set_ylim(axes[2].get_ylim())

    ## axis label
    axes[2].xaxis.set_major_locator(FixedLocator(np.arange(1, 2.51, 0.5)))  # FixedLocator
    axes[2].xaxis.set_minor_locator(FixedLocator(np.arange(1, 2.51, 0.1)))  # FixedLocator
    axes[2].yaxis.set_major_locator(FixedLocator(np.arange(-50, 40.01, 10)))  # FixedLocator
    axes[2].yaxis.set_minor_locator(FixedLocator(np.arange(-50, 40.01, 2)))  # FixedLocator

    # --------------------------------
    ## add subplot label
    labels = 'AB'.lower()
    for i, label in enumerate(labels):
        if i == 0:
            xloc = -0.20;
            yloc = 1.05
            ax = axes[0]

        if i == 1:
            xloc = -0.20;
            yloc = 1.05
            ax = axes[2]

        ax.text(xloc, yloc, label, transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


    ## save figure
    fig.savefig(f'manuscript figures/pdf/{figname}.pdf')
    fig.savefig(f'manuscript figures/png/{figname}.png')
    
    fig.clf()
    plt.cla()


    print("done")
    return df_2Dgrid_mean



if __name__ == "__main__":
    import os
    import argparse


    if os.getcwd().find("manuscript_figure_script") == -1:
        parser = argparse.ArgumentParser()
        parser.add_argument('--figname', type=str, default="figure3")
        parser.add_argument('--llk_2dgrid', type=str, required=True)
        parser.add_argument('--phydyn_meta', type=str, required=False, default=None)
        parser.add_argument('--phydyn_log', type=str, required=False, default=None)

        args = parser.parse_args()

    else:
        os.chdir(os.getcwd().replace("manuscript_figure_script", ""))
        parser = argparse.ArgumentParser()
        args = parser.parse_args()
        args.llk_2dgrid  = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_gridsearch_R0,timestart"
        args.phydyn_meta = None
        args.phydyn_log  = None
        args.figname     = "figure3"



    print(args)

    hey = run(args.llk_2dgrid, args.phydyn_meta, args.phydyn_log, args.figname)
