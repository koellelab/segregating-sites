import os
import pickle
import datetime
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator

from fig2_plot import get_MLE_and_CI, read_and_get_mean_llk
from fig3_plot import plot_joint_estimation_95CI
# ---------------------------
def run(args):
    txt = open(f'manuscript figures/txt/{args.figname}.txt', 'w')
    # ---------------------------
    plt.rcParams.update({'font.size': 10})
    fig = plt.figure(figsize=[10, 2.85], dpi=300)

    gs0 = fig.add_gridspec(1, 3)
    axes = [fig.add_subplot(gs0[0]),
            fig.add_subplot(gs0[1]),
            fig.add_subplot(gs0[2])]
    # ---------------------------
    for i, dir_llk_2Dgrid in enumerate(args.llk_2dgrid):
        df_2Dgrid_mean, _ = read_and_get_mean_llk(dir_llk_2Dgrid, args.params, args.n_iter_per_cell, return_original = True)


        print(df_2Dgrid_mean.sort_values('mean_logL', ascending=False).iloc[:10].to_string())
        df_2Dgrid_mean[df_2Dgrid_mean.sim_no == args.n_iter_per_cell].to_csv(dir_llk_2Dgrid + "/meanllk.tsv", sep="\t", index=False)

        df_2Dgrid_mean.mean_logL = np.where(~np.isnan(df_2Dgrid_mean.mean_logL), df_2Dgrid_mean.mean_logL, 1e5)
        df_2Dgrid_mean.mean_logL = np.where(~np.isinf(df_2Dgrid_mean.mean_logL), df_2Dgrid_mean.mean_logL, 1e5)
        df_2Dgrid_mean.mean_logL = np.where(df_2Dgrid_mean.sim_no >= args.n_iter_per_cell, df_2Dgrid_mean.mean_logL, np.nan)

        plot_joint_estimation_95CI(axes[i], fig, df_2Dgrid_mean, ['R0', 'log10eta'], [np.nan, np.nan], gamma=55, majortick=5, majortick_y=5, df = 2)

        ## axis
        axes[i].set_yticklabels(["%.1f" % (np.unique(df_2Dgrid_mean.log10eta)[j])
                                    for j in range(len(np.unique(df_2Dgrid_mean.log10eta))) if j % 5 == 0])
        axes[i].set_xlabel(r"$R_0$")
        axes[i].set_ylabel(r"log$_{10}\eta$")

        # print MLE and CI to txt
        df_2Dgrid_mean.mean_logL = np.where(df_2Dgrid_mean.mean_logL != 1e5, df_2Dgrid_mean.mean_logL, np.nan)
        mle_and_CI = get_MLE_and_CI(df_2Dgrid_mean, 'R0', MA=False, df = 2)
        print("--- R0 ----\n", file=txt)
        print(mle_and_CI[0], "\n", file=txt)
        print(mle_and_CI[1], file=txt)
        print(mle_and_CI[2], "\n", file=txt)
        print(mle_and_CI[1], "-", mle_and_CI[2])

        mle_and_CI = get_MLE_and_CI(df_2Dgrid_mean, 'log10eta', MA=False, df = 2)
        print("--- eta ----\n", file=txt)
        print(mle_and_CI[0], "\n", file=txt)
        print(mle_and_CI[1], file=txt)
        print(mle_and_CI[2], "\n", file=txt)
        print(f"{mle_and_CI[0]['log10eta']:.3f}, {mle_and_CI[1]:.3f} - {mle_and_CI[2]:.3f}", file=txt)
        print ("===========================\n", file=txt)















    txt.close()
    # ---------------------
    ## add subplot label
    label = "ABC".lower()
    for i in range(len(label)):
        xloc = -0.27;
        yloc = 1.1

        axes[i].text(xloc, yloc, label[i], transform=axes[i].transAxes, fontsize=14, fontweight='bold', va='top',
                     ha='left')

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
        parser.add_argument('--figname', type=str, default="figure8")
        parser.add_argument('--n_iter_per_cell', type=int, default=20)
        parser.add_argument('--params', type=str, nargs=2,required=True)
        parser.add_argument('--llk_2dgrid', type=str, nargs=3, required=True)

        args = parser.parse_args()
        print(args)


    else:
        os.chdir(os.getcwd().replace("manuscript_figure_script", ""))
        parser = argparse.ArgumentParser()
        args = parser.parse_args()

        print(os.getcwd())
        args.figname = "figure8"
        args.n_iter_per_cell = 10
        args.llk_2dgrid = ["empirical_data/france_data/france_multiple_te=191224/gridsearch_log10eta,R0",
                           "empirical_data/france_data/france_multiple_te=200101/gridsearch_log10eta,R0",
                           "empirical_data/france_data/france_multiple_te=200108/gridsearch_log10eta,R0"]
        args.params = ["R0", "log10eta"]
#"/Users/yeongseon/Dropbox/PhD_Emory_university/1_projects/Project_seg_site/segregating-sites/empirical_data/france_data/france_multiple_te=191224/gridsearch_log10eta,R0"



    run(args)





