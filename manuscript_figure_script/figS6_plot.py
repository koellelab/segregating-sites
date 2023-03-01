import os
import datetime
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator

from fig1_plot import plot_data, plot_statevar
from figS4_plot import n_sequence_over_time, n_segregating_over_time, read_and_get_mean_llk, get_MLE_and_CI, plot_MLE_and_CI
# ---------------------------
def run (args):
    dirnames = list(zip(args.simul_data, args.llk_rand))
    txt = open(f'manuscript figures/txt/{args.figname}.txt', 'w')
    # ---------------------------
    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(1 * 3, 5, figsize=[8 * 1.4, 2.25 * 3], dpi=300)
    # ---------------------------
    for i, dirname in enumerate(dirnames):
        dir_simulated_dataset, dir_likelihood_randsearch = dirname

        ## (A) number of samples per window
        ## (B) number of segregating sites per window (both data and reconstructed)
        df_simulated_dataset = pd.read_csv(dir_simulated_dataset, sep="\t")
        n_sequence_over_time   (axes[0, i], data=df_simulated_dataset)
        n_segregating_over_time(axes[1, i], data=df_simulated_dataset)

        ## (C) Likelihood at t0 = 0
        df_randsearch, df = read_and_get_mean_llk(dir_likelihood_randsearch, ['R0'], args.n_iter_per_cell, return_original=True)  # 19np.inf)
        df_randsearch.to_csv(dir_likelihood_randsearch + "/meanllk.tsv", sep="\t", index=False)

        axes[2, i].scatter(df['R0'], df['logL'], edgecolors="grey", s=6, facecolors='none', alpha=0.3)
        axes[2, i].plot(df_randsearch['R0'], df_randsearch['mean_logL'], linewidth=1, c="k")

        mle_and_CI = get_MLE_and_CI(df_randsearch, 'R0', MA=False)
        plot_MLE_and_CI(axes[2, i], mle_and_CI, 'R0')
        axes[2, i].axvline(x=args.true_param, linewidth=1, ls='--')  ## true value

        # ---------------------------
        ## set axis
        axes[2, i].set_xlabel("Days post index case")
        axes[2, i].set_xlim([0.95, 2.55])
        axes[2 ,i].xaxis.set_major_locator(FixedLocator(np.arange(1.0, 2.5, 0.5).tolist() + [2.5]))
        axes[2 ,i].xaxis.set_minor_locator(FixedLocator(np.arange(1.0, 2.5, 0.1).tolist() + [2.5]))


        # print MLE and CI to txt
        print(f'run at {datetime.datetime.now().strftime("%y%m%d %H:%M:%S")}\n')
        print(mle_and_CI[0], "\n", file=txt)
        print(mle_and_CI[1], file=txt)
        print(mle_and_CI[2], "\n", file=txt)
        print(f"{mle_and_CI[1]:.4f} - {mle_and_CI[2]:.4f}", file=txt)

    # ticks and labels
    axes[0, 0].set_ylabel("Number of\nsequences")
    axes[1, 0].set_ylabel("Number of\nsegregating sites")
    axes[2, 0].set_ylabel('log-likelihood')

    ## add subplot label
    labels = np.array([x for x in 'ABCDEFGHIJKLMNO'.lower()]).reshape(5, 3)
    labels = ["".join(labels[:,i]) for i in range(labels.shape[1])]
    #labels = ['ABCDE', 'FGHIJ', 'KLMNO']
    for i, label_row in enumerate(labels):
        for j, l in enumerate(label_row):
            if j == 0:
                x_loc = -0.4
                axes[i][j].text(-0.35, 1.1, l, transform=axes[i][j].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')
            else:
                x_loc = -0.3
                axes[i][j].text(-0.23, 1.1, l, transform=axes[i][j].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')
            continue

    txt.close()

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
        parser.add_argument('--figname', type=str, default="figure3")
        parser.add_argument('--simul_data', type=str, required=True, nargs='+')
        parser.add_argument('--llk_rand', type=str, required=True, nargs='+')
        parser.add_argument('--true_param', type=float, required=True)
        parser.add_argument('--xlim', type=float, nargs=2)
        parser.add_argument('--ylim', type=float, nargs=2)
        parser.add_argument('--n_iter_per_cell', type=int)

        args = parser.parse_args()


    else:
        os.chdir(os.getcwd().replace("manuscript_figure_script", ""))
        parser = argparse.ArgumentParser()
        args = parser.parse_args()


    if len(args.simul_data) != len(args.llk_rand):
        raise Exception()


    print(args)
    hey = run(args)



