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
from datetime import datetime
from datetime import timedelta

def convert_to_matlab_date (datestr, sep = "/"):
    return datetime.strptime(datestr, sep.join(['%Y', '%m', '%d'])).toordinal() + 366

def convert_to_caldate (matlab_datenum, sep = "/"):
    return (datetime(1, 1, 1) + timedelta(matlab_datenum - 367)).strftime(sep.join(['%Y', '\n%m', '%d']))



def run(args):
    dir_llk_2Dgrid = args.llk_2dgrid
    # ---------------------------
    plt.rcParams.update({'font.size': 10})
    fig = plt.figure(figsize=[5, 4], dpi = 300)
    gs0 = fig.add_gridspec(1, 1)
    axes = [fig.add_subplot(gs0[0])]
    # ---------------------------
    ## panel A
    df = 2
    df_2Dgrid_mean = read_and_get_mean_llk(dir_llk_2Dgrid, args.params, args.n_iter_per_cell)
    df_2Dgrid_mean.to_csv(dir_llk_2Dgrid + "/meanllk.tsv", sep="\t", index=False)

    df_2Dgrid_mean.mean_logL = np.where(~np.isinf(df_2Dgrid_mean.mean_logL), df_2Dgrid_mean.mean_logL, 1e5)
    df_2Dgrid_mean.mean_logL = np.where(~np.isnan(df_2Dgrid_mean.mean_logL), df_2Dgrid_mean.mean_logL, 1e5)
    df_2Dgrid_mean.mean_logL = np.where(df_2Dgrid_mean.sim_no == args.n_iter_per_cell, df_2Dgrid_mean.mean_logL, np.nan)

    plot_joint_estimation_95CI(axes[0], fig, df_2Dgrid_mean, ["R0", "timestart"], [np.nan, np.nan], gamma=20,  df = 2, majortick=5)


    hey = df_2Dgrid_mean.pivot("timestart", "R0", values='sim_no')



    # ---------------------------
    ## label figure
    for i in range(1):
        axes[i].set_xlabel(r'$R_0$', fontsize=14)
        axes[i].set_ylabel(r'$t_0$', fontsize=14)

        ## convert matlab date to caldate
        yticklabel = axes[i].get_yticklabels()
        yticklabel = [convert_to_caldate(float(x.get_text())).replace("2018/", "") for x in yticklabel]
        axes[i].set_yticklabels(yticklabel)

    # print MLE and CI to txt; R0
    df_2Dgrid_mean.mean_logL = np.where(df_2Dgrid_mean.mean_logL != 1e5, df_2Dgrid_mean.mean_logL, np.nan)

    txt = open(f'manuscript figures/txt/{args.figname}.txt', 'w')
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
        parser.add_argument('--figname', type=str, default="figure7")
        parser.add_argument('--llk_2dgrid', type=str, required=True)
        parser.add_argument('--params', type=str, nargs=2,required=True)
        parser.add_argument('--n_iter_per_cell', type=int, default=20)


        args = parser.parse_args()


    else:
        os.chdir(os.getcwd().replace("manuscript_figure_script", ""))
        parser = argparse.ArgumentParser()
        args = parser.parse_args()

        args.figname = "figure7"
        args.params = ["timestart", "R0"]
        args.llk_2dgrid  = "empirical_data/france_data/france_single/gridsearch_timestart,R0"
        args.n_iter_per_cell = 10





    #print(args)
    run(args)

