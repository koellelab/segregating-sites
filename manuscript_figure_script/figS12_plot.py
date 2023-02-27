import os
import pickle
import datetime
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator

from fig7_plot import convert_to_matlab_date, convert_to_caldate
# ---------------------------

def run(args):
    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=[8, 2.85], dpi=300)
    # -----------------------------
    df = [pd.read_csv(x, sep ="\t") for x in args.dir_data]

    for i in range(len(df)):
        df[i]['window_end_cal'] = [convert_to_caldate(x).replace("20d20/\n", "") for x in df[i]['window_end']]
        axes[0].plot(df[i]['window_end_cal'], df[i]['s'])
        axes[0].scatter(df[i]['window_end_cal'], df[i]['s'], s=15)

        if i != 2:
            axes[1].plot(df[i]['window_end_cal'], df[i]['s'], linewidth=1)
            axes[1].scatter(df[i]['window_end_cal'], df[i]['s'], s=15)


    axes[0].set_ylim([-20, 2000])
    axes[1].set_ylim([-2, 200])

    axes[0].set_xlabel("Days post index case")
    axes[0].set_ylabel("Number of segregating sites")
    axes[1].set_xlabel("Days post index case")
    axes[1].set_ylabel("Number of segregating sites")

    ## save figure
    plt.tight_layout(pad=0.5)
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
        parser.add_argument('--figname', type=str, default="figureS12")
        parser.add_argument('--dir_data', nargs=3, type=str)

        args = parser.parse_args()

    else:
        os.chdir(os.getcwd().replace("manuscript_figure_script", ""))
        parser = argparse.ArgumentParser()
        args = parser.parse_args()

    print(args)
    run(args)
