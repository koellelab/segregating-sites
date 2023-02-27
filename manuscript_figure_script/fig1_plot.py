import os
import datetime
import numpy as np
import pandas as pd
    
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator
# ---------------------------
def plot_data(ax, label, c, filename, n_trial=30, time_shift=0,  ls="-", alpha=0.2, linewidth = None):
    for trial in range(n_trial):
        print(filename.replace("##", str(trial)))
        data = pd.read_csv(filename.replace("##", str(trial)), sep="\t")
        data = data.loc[~np.isnan(data['s'])]
        if linewidth is None:
            ax.plot(data['window_end']+time_shift, data['s'], alpha=alpha, label=label, c=c, ls=ls)
        else:
            ax.plot(data['window_end']+time_shift, data['s'], alpha=alpha, label=label, c=c, ls=ls, linewidth=linewidth)
        label = None
    return

def plot_statevar(ax, label, c, filename, time_shift=0,  ls="-", alpha=0.7, names = ["t", "S", "E", "I", "R"]):
    statevar = pd.read_csv(filename, sep="\t", comment="#", names=names)
    if "I" in names:
        ax.plot(statevar['t']+time_shift, statevar['I'], alpha=alpha, label=label, c=c, ls=ls)
    else:
        print(statevar)
        ax.plot(statevar['t'] + time_shift, statevar['I_l'] + statevar['I_h'], alpha=alpha, label=label, c=c, ls=ls)



def run(baseline_statevar,baseline_dense,baseline_sparse,baseline_full,
        highR0, highR0_dense, transhet, transhet_dense,
        reduced1p1, reduced1p1_dense, reduced0p75, reduced0p75_dense, figname):

    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(2, 4, figsize=[13.5, 5], dpi=300)

    # ---------------------------
    ## Panel A, B: dense vs. sparse vs. full
    plot_statevar(axes[0, 0], 'dense',   "k",       baseline_statevar)
    plot_data    (axes[1, 0], 'dense',   "k",       baseline_dense)
    plot_data    (axes[1, 0], 'sparse',  "#9e9e9e", baseline_sparse)
    plot_data    (axes[1, 0], 'full',    "k"      , baseline_full, ls = "--", n_trial = 1, alpha = 1)

    ## inset - same data but zoomed out
    left, bottom, width, height = [0.07, 0.33, 0.075, 0.12]
    inset = fig.add_axes([left, bottom, width, height])

    plot_data    (inset, 'dense',   "k",       baseline_dense)
    plot_data    (inset, 'sparse',  "#9e9e9e", baseline_sparse)
    plot_data    (inset, 'full',    "k"      , baseline_full, ls = "--", n_trial = 1, alpha = 1)

    inset.patch.set_alpha(0.8)
    inset.set_xlim([-25, 205])
    inset.xaxis.set_ticklabels([])

    inset.yaxis.set_major_locator(FixedLocator([]))
    inset.set_ylim([-10, 2750])
    inset.yaxis.set_major_locator(FixedLocator([200, 1500, 2500]))

    inset.tick_params(axis = 'y', direction = 'in', pad = -67, labelsize = 6, labelright = True, labelleft = False, right = False, left = True)


    # ---------------------------
    ## Panel C, D: high R0 vs. low R0 - todo
    plot_statevar(axes[0, 1], r'R$_{0}$ = 1.6', "k",       baseline_statevar)
    plot_data    (axes[1, 1], r'R$_{0}$ = 1.6', "k",       baseline_dense)

    plot_statevar(axes[0, 1], r'R$_{0}$ = 2.0', "#318fdb", highR0)
    plot_data    (axes[1, 1], r'R$_{0}$ = 2.0', "#318fdb", highR0_dense)



    # ---------------------------
    # Panel E, F: transmission heterogeneity
    statevar_1 = pd.read_csv(baseline_statevar, sep="\t", comment="#", names=["t", "S", "E",  "I", "R"])
    statevar_2 = pd.read_csv(transhet, sep="\t", comment="#",names=["t", "S", "E", "I_l", "I_h", "R"])
    time_shift =  statevar_1['t'][np.argmax(statevar_1['I'])] - \
                  statevar_2['t'][np.argmax(statevar_2['I_l'] + statevar_2['I_h'])]

    plot_statevar(axes[0, 2], 'No TH',       "k",       baseline_statevar)
    plot_data    (axes[1, 2], 'No TH',       "k",       baseline_dense)

    plot_statevar(axes[0, 2], 'TH'          , "#80b564", transhet,       ls = "-.", names = ["t", "S", "E", "I_l", "I_h", "R"])
    plot_statevar(axes[0, 2], 'TH (shifted)', "#80b564", transhet,       time_shift= time_shift, names = ["t", "S", "E", "I_l", "I_h", "R"])
    plot_data    (axes[1, 2], 'TH (shifted)', "#80b564", transhet_dense, time_shift= time_shift)



    # ---------------------------
    # ## Panel G, H: reduced R0
    plot_statevar(axes[0, 3], r'R$_{0}$ unchanged',       "k",       baseline_statevar)
    plot_data    (axes[1, 3], r'R$_{0}$ unchanged',       "k",       baseline_dense)

    plot_statevar(axes[0, 3], r'R$_{0}$ reduced to 1.1' , "#ff5f43", reduced1p1)
    plot_data    (axes[1, 3], r'R$_{0}$ reduced to 1.1' , "#ff5f43", reduced1p1_dense)
    plot_statevar(axes[0, 3], r'R$_{0}$ reduced to 0.75', "#fcb708", reduced0p75)
    plot_data    (axes[1, 3], r'R$_{0}$ reduced to 0.75', "#fcb708", reduced0p75_dense)
    #

    # ---------------------------
    # set axis limit
    for i in range(4):
        axes[0, i].set_xlim([0, 200]); axes[0, i].set_ylim([0, 10000]);
        axes[1, i].set_xlim([0, 200]); axes[1, i].set_ylim([0, 225]);
        axes[0, i].ticklabel_format(axis='y', style="sci", scilimits=(3, 4))

        if i != 0:
            axes[0, i].legend(loc="upper right")
        axes[1, i].set_xlabel('Days post index case')

    axes[0, 0].set_ylabel("Number of\ninfected individuals")
    axes[1, 0].set_ylabel("Number of\nsegregating sites")


    ## add subplot label
    labels = ['ACEG'.lower(), 'BDFH'.lower()]
    for i, label_row in enumerate(labels):
        for j, l in enumerate(label_row):
            if j == 0:
                x_loc = -0.4
                axes[i][j].text(-0.27, 1.1, l, transform=axes[i][j].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')
            else:
                x_loc = -0.3
                axes[i][j].text(-0.23, 1.1, l, transform=axes[i][j].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')
            continue


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
        parser.add_argument('--figname', type=str, default="figure3")

        parser.add_argument('--baseline_statevar',  type=str, required=True)
        parser.add_argument('--baseline_dense',     type=str, required=True)
        parser.add_argument('--baseline_sparse',    type=str, required=True)
        parser.add_argument('--baseline_full',      type=str, required=True)

        parser.add_argument('--highR0',         type=str, required=True)
        parser.add_argument('--highR0_dense',   type=str, required=True)

        parser.add_argument('--transhet',       type=str, required=True)
        parser.add_argument('--transhet_dense', type=str, required=True)

        parser.add_argument('--reduced1p1',         type=str, required=True)
        parser.add_argument('--reduced1p1_dense',   type=str, required=True)
        parser.add_argument('--reduced0p75',        type=str, required=True)
        parser.add_argument('--reduced0p75_dense',  type=str, required=True)

        args = parser.parse_args()

    else:
        os.chdir(os.getcwd().replace("manuscript_figure_script", ""))
        parser = argparse.ArgumentParser()
        args = parser.parse_args()
        args.figname = "figure1"

        args.baseline_statevar = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed20230201_statevar.tsv"
        args.baseline_dense    = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed20230201_unif40_win4/seed230207_segsites.tsv"
        args.baseline_sparse   = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed20230201_unif20_win4/seed230207_segsites.tsv"
        args.baseline_full     = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed20230201_full_win4/seed230201_segsites.tsv"

        args.highR0            = "simulated_data"
        args.highR0_dense      = "simulated_data"

        args.transhet          = "simul_dat"
        args.transhet_dense    = "simul_dat"

        args.reduced1p1        = "simul_dat"
        args.reduced1p1_dense  = "simul_dat"
        args.reduced0p75       = "simul_dat"
        args.reduced0p75_dense = "simul_dat"


    print(args)
    run(args.baseline_statevar,
        args.baseline_dense,
        args.baseline_sparse,
        args.baseline_full,

        args.highR0,
        args.highR0_dense,
        args.transhet,
        args.transhet_dense,

        args.reduced1p1,
        args.reduced1p1_dense,
        args.reduced0p75,
        args.reduced0p75_dense,

        args.figname)


