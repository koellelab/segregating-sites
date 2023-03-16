import os
import pickle
import datetime
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator
# ---------------------------

def run(args):
    dir_high    = args.dir_high
    dir_default = args.dir_default
    # ---------------------------
    mu_high = pd.read_csv(dir_high, sep="\t")
    default = np.load(dir_default, allow_pickle=True)['statevar']


    plt.plot(default[:, 0], default[:, 2], label='simulation for 5A')
    plt.plot(mu_high.iloc[:, 0], mu_high.I, label='simulation for 5D')
    plt.xlabel('Days post index case')
    plt.ylabel('Number of infected individuals')

    plt.xlim([0, 70])
    plt.ylim([0, 100])

    plt.legend()

    ## save figure
    # plt.tight_layout(pad=0.5)
    plt.savefig(f'manuscript figures/pdf/{args.figname}.pdf')
    plt.savefig(f'manuscript figures/png/{args.figname}.png')

    # fig.clf()
    plt.cla()

    return



if __name__ == "__main__":
    import os
    import argparse


    if os.getcwd().find("manuscript_figure_script") == -1:
        parser = argparse.ArgumentParser()
        parser.add_argument('--figname', type=str, default="figureS8")
        parser.add_argument('--dir_high', type=str)
        parser.add_argument('--dir_default', type=str)


        args = parser.parse_args()

    else:
        os.chdir(os.getcwd().replace("manuscript_figure_script", ""))
        parser = argparse.ArgumentParser()
        args = parser.parse_args()

    print(args)
    run(args)
