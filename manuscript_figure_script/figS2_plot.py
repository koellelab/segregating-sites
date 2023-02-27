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
    color = prop_cycle.by_key()['color']
    # ---------------------------
    statevar_1 = pd.read_csv(args.statevar_080, sep="\t", comment="#", names=["t", "S", "E", "I", "R"])
    plot_statevar(axes[0], r'$p_{H}$=0.8',       "k", args.statevar_080, alpha = 0.7)
    plot_data    (axes[1], r'$p_{H}$=0.8',       "k", args.segsite_080, linewidth = 1)


    fnames = [[args.statevar_006, args.segsite_006],
              [args.statevar_030, args.segsite_030],
              [args.statevar_015, args.segsite_015]]

    for i, p_H in enumerate([0.06, 0.15, 0.3]):
        statevar_2 = pd.read_csv(fnames[i][0], sep="\t", comment="#",names=["t", "S", "E", "I_l", "I_h", "R"])
        time_shift =    statevar_1['t'][np.argmax(statevar_1['I'])] - \
                        statevar_2['t'][np.argmax(statevar_2['I_l'] + statevar_2['I_h'])]

        plot_statevar(axes[0], r'$p_{H}$='+str(p_H), color[i], fnames[i][0], time_shift= time_shift, names=["t", "S", "E", "I_l", "I_h", "R"], alpha = 0.7)
        plot_data    (axes[1], r'$p_{H}$='+str(p_H), color[i], fnames[i][1], time_shift= time_shift, linewidth = 1)



    # ---------------------------
    ## xticks
    for i in [0, 1]:
        axes[i].xaxis.set_major_locator(FixedLocator(np.arange(00, 200.1, 24)))  # FixedLocator
        axes[i].xaxis.set_minor_locator(FixedLocator(np.arange(00, 200.1, 8)))  # FixedLocator
        axes[i].set_xlim([0, 200])
        axes[i].grid(which="both")

    # yticks
    axes[0].set_ylim([0, 10000]);
    axes[0].ticklabel_format(axis='y', style="sci", scilimits=(3, 4))
    axes[1].set_ylim([0, 225]);

    ## legend
    handles, labels = plt.gca().get_legend_handles_labels()
    handles, labels = axes[0].get_legend_handles_labels()
    print(axes[0].get_legend_handles_labels())
    order = [0, 3, 2, 1]
    axes[0].legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc="upper right")
    # axes[0].legend()

    # axis label
    axes[0].set_ylabel("Number of infected individuals")
    axes[1].set_ylabel("Number of segregating sites")
    axes[1].set_xlabel('Days post index case')

    ## add subplot label
    labels = 'AB'.lower()
    for i, label in enumerate(labels):
        axes[i].text(-0.22, 1.1, label, transform=axes[i].transAxes, fontsize=16, fontweight='bold', va='top',
                     ha='left')

    ## save figure
    fig.tight_layout(pad=0.5)
    fig.savefig(f'manuscript figures/pdf/{args.figname}.pdf')
    fig.savefig(f'manuscript figures/png/{args.figname}.png')

    fig.clf()
    plt.cla()





if __name__ == "__main__":
    import os
    import argparse


    if os.getcwd().find("manuscript_figure_script") == -1:
        parser = argparse.ArgumentParser()
        parser.add_argument('--figname', type=str, default="figureS2")
        parser.add_argument('--statevar_080', type=str)
        parser.add_argument('--segsite_080' , type=str)
        parser.add_argument('--statevar_006', type=str)
        parser.add_argument('--segsite_006' , type=str)
        parser.add_argument('--statevar_015', type=str)
        parser.add_argument('--segsite_015' , type=str)
        parser.add_argument('--statevar_030', type=str)
        parser.add_argument('--segsite_030' , type=str)

        args = parser.parse_args()

    else:
        os.chdir(os.getcwd().replace("manuscript_figure_script", ""))
        parser = argparse.ArgumentParser()
        args = parser.parse_args()
        args.figname = "figureS2"

        args.statevar_080 = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed20230201_statevar.tsv"
        args.segsite_080  = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed20230201_unif40_win4/seed230201_trial##_segsites.tsv"
        args.statevar_006 = "simulated_data/heteroSEIR_R0=16e-1_mu2e-1/seed20230201_statevar.tsv"
        args.segsite_006  = "simulated_data/heteroSEIR_R0=16e-1_mu2e-1/seed20230201_unif40_win4/seed230201_trial##_segsites.tsv"
        args.statevar_015 = "simulated_data/heteroSEIR_R0=16e-1_mu2e-1_ph=15e-2/seed20230201_statevar.tsv"
        args.segsite_015  = "simulated_data/heteroSEIR_R0=16e-1_mu2e-1_ph=15e-2/seed20230201_unif40_win4/seed230201_trial##_segsites.tsv"
        args.statevar_030 = "simulated_data/heteroSEIR_R0=16e-1_mu2e-1_ph=3e-1/seed20230201_statevar.tsv"
        args.segsite_030  = "simulated_data/heteroSEIR_R0=16e-1_mu2e-1_ph=3e-1/seed20230201_unif40_win4/seed230201_trial##_segsites.tsv"


    print(args)
    run(args)

