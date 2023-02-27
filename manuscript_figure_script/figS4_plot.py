import os
import pickle
import datetime
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator
# ---------------------------
def read_and_get_mean_llk (dir_llk_grid2D, params, n_per_cell, return_original=False):
    fnames_likelihood_randsearch = [x for x in os.listdir(dir_llk_grid2D) if
                                        (x.find(".pkl") > -1) and
                                        (not x.find('window.pkl') > -1)]
    ## read files
    df_2Dgrid = pd.DataFrame({})
    for f in fnames_likelihood_randsearch:
        #print(f)
        tmp = pd.DataFrame(pickle.load(open(dir_llk_grid2D + "/" + f, "rb")))
        tmp.columns = ['sim_no'] + params + ['logL']
#        if len(tmp.columns) == 4:
#            tmp.columns = ['sim_no'] + params + ['logL']
#        elif len(tmp.columns) == 3:
#            tmp.columns = params + ['logL']
#            tmp['sim_no'] = 0

        tmp.R0 = np.round(tmp.R0 * 1e6) / 1e6
        tmp = tmp[tmp.sum(axis=1) != 0]

        df_2Dgrid = pd.concat([df_2Dgrid, tmp])

    print("")

    if len(params) > 2 and params[1] == "eta":
        df_2Dgrid['log10eta'] = np.log10(df_2Dgrid.eta)
        params = ['R0', 'log10eta']

    ## assign id number to filter out excess iterations
    df_group = df_2Dgrid.groupby(params, as_index=False)
    df_2Dgrid['id'] = df_group[params[0]].rank(method='first') - 1
    df_2Dgrid = df_2Dgrid.sort_values(params+['id'], ascending=False)
    df_2Dgrid = df_2Dgrid[df_2Dgrid.id < n_per_cell]

    ## group by parameter sets
    df_group_mean = df_2Dgrid.groupby(params, as_index=False)\
        .agg({'logL': 'mean', 'sim_no': 'size'}) \
        .rename(columns={'logL': 'mean_logL'})

    print (f">> {np.unique(df_group_mean.sim_no)} iterations per cell")

    if return_original:
        return df_group_mean, df_2Dgrid
    else:
        return df_group_mean


def n_sequence_over_time(ax, data):

    data = data.loc[~np.isnan(data['n'])]
    ax.plot(data['window_end'], data['n'], c="k")
    ax.scatter(data['window_end'], data['n'], c="k", s=15)

    #ax.set_xlabel("Days post index case")
    #ax.set_ylabel("Number of sequences")

    return

def n_segregating_over_time(ax, data):

    data = data.loc[~np.isnan(data['s'])]
    print(data)
    print("\n\n\n\n")

    ax.plot(data['window_end'], data['s'], c="k")
    ax.scatter(data['window_end'], data['s'], c="k", s=15)

    #ax.set_xlabel("Days post index case")
    #ax.set_ylabel("Number of segregating sites")

    return

def sampling_prop_over_time(ax, data):

    data = data.loc[~np.isnan(data['s'])]
    print(data)
    print("\n\n\n\n")

    ax.plot   (data['window_end'], data['n'] / data['n_recovered_in_win'], c="k", linestyle = "--", linewidth = 0.7)
    ax.scatter(data['window_end'], data['n'] / data['n_recovered_in_win'], c="k", s=8, marker = "s")

    #ax.set_xlabel("Days post index case")
    #ax.set_ylabel("Number of segregating sites")

    return

def get_MLE_and_CI(df_randsearch, param, MA = True, df = 1):
    
    if df == 1:
        thres = 1.92
    if df == 2:
        thres = 2.9957325
    
    if MA:
        max_MA = df_randsearch.iloc[np.argmax(df_randsearch['MA'])]

        CI = df_randsearch[df_randsearch.MA > max_MA['MA'] -thres]
        CI_logL_lower = CI[param].min()  # CI_logL.iloc[np.argmin(CI_logL['R0'])]
        CI_logL_upper = CI[param].max()

    else:
        max_MA = df_randsearch.iloc[np.argmax(df_randsearch['mean_logL'])]

        CI = df_randsearch[df_randsearch.mean_logL > max_MA['mean_logL'] - thres]
        CI_logL_lower = CI[param].min()  # CI_logL.iloc[np.argmin(CI_logL['R0'])]
        CI_logL_upper = CI[param].max()


    return max_MA, CI_logL_lower, CI_logL_upper, CI

def plot_MLE_and_CI(ax, MLE_and_CI, param, df = 1):

    if df == 1:
        thres = 1.92
    elif df == 2:
        thres = 2.9957325


    max_MA, CI_logL_lower, CI_logL_upper, _ = MLE_and_CI

    ax.axvline(x=max_MA[param], linewidth=1, ls='-.', color="firebrick")
    ax.axvspan(CI_logL_lower, CI_logL_upper, alpha=0.1, color='red')
    ax.axhspan(max_MA['mean_logL'] - thres, max_MA['mean_logL'], alpha=0.1, color='steelblue')

    return




def run(dir_simulated_dataset, dir_likelihood_randsearch, true_param, xlim, ylim, figname):

    # ---------------------------
    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(1, 3, figsize=[8, 2.85], dpi=300)
    # ---------------------------
    ## (A) number of samples per window
    ## (B) number of segregating sites per window (both data and reconstructed)
    df_simulated_dataset = pd.read_csv(dir_simulated_dataset, sep="\t")
    n_sequence_over_time   (axes[0], data=df_simulated_dataset)
    n_segregating_over_time(axes[1], data=df_simulated_dataset)

    ## (C) Likelihood at t0 = 0
    n_iter_per_cell=20
    df_randsearch, df = read_and_get_mean_llk(dir_likelihood_randsearch, ['R0'], n_iter_per_cell, return_original=True)  # 19np.inf)
    axes[2].scatter(df['R0'], df['logL'], edgecolors="grey", s=6, facecolors='none', alpha=0.3)
    axes[2].plot   (df_randsearch['R0'], df_randsearch['mean_logL'], linewidth=1, c="k")

    mle_and_CI = get_MLE_and_CI(df_randsearch, 'R0', MA = False)
    plot_MLE_and_CI(axes[2], mle_and_CI, 'R0')
    axes[2].axvline(x=true_param, linewidth=1, ls='--')        ## true value


    # ---------------------------
    ## set axis
    axes[0].set_xlabel("Days post index case")
    axes[0].set_ylabel("Number of\nsequences")
    axes[1].set_xlabel("Days post index case")
    axes[1].set_ylabel("Number of\nsegregating sites")
    axes[2].set_xlabel(r'$R_0$')
    axes[2].set_ylabel('log-likelihood')

    axes[0].set_ylim([-0.3, int(df_simulated_dataset['n'].max() * 1.1)])

    if np.unique(df_simulated_dataset['n']).size == 1:
        axes[0].set_ylim([-0.3, int(df_simulated_dataset['n'].max() * 1.7)])
        axes[0].set_yticks([0, df_simulated_dataset['n'].max()])
    
    if ylim:
        axes[2].set_ylim([ylim[0], ylim[1]])

    if xlim:
        axes[2].set_xlim([xlim[0], xlim[1]])
        axes[2].xaxis.set_major_locator(FixedLocator(np.arange(xlim[0], xlim[1], 0.5).tolist() + [xlim[1]]))
        axes[2].xaxis.set_minor_locator(FixedLocator(np.arange(xlim[0], xlim[1], 0.1).tolist() + [xlim[1]]))

    
    ## add subplot label
    labels = 'ABC'.lower()
    for i, l in enumerate(labels):
        if i == 0:
            axes[i].text(-0.42, 1.1, l, transform=axes[i].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')
        else:
            axes[i].text(-0.45, 1.1, l, transform=axes[i].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


    # print MLE and CI to txt
    txt = open(f'manuscript figures/txt/{figname}.txt', 'w')
    print(f'run at {datetime.datetime.now().strftime("%y%m%d %H:%M:%S")}\n')

    #print(f">> using {n1} out of  {n2} data points (ranges from R0 = {df_randsearch['R0'].min():.4f} to {df_randsearch['R0'].max():.4f})", file=txt)
    #print(f">> window size = {n3}\n", file = txt)
    print(mle_and_CI[0], "\n",  file = txt)
    print(mle_and_CI[1], file = txt)
    print(mle_and_CI[2], "\n", file = txt)
    print(f"{mle_and_CI[1]:.4f} - {mle_and_CI[2]:.4f}", file = txt)

    txt.close()

    ## save figure
    fig.suptitle(figname)
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
        parser.add_argument('--simul_data', type=str, required=True)
        parser.add_argument('--llk_rand', type=str, required=True)
        parser.add_argument('--true_param', type=float, required=True)
        parser.add_argument('--xlim', type=float, nargs=2)
        parser.add_argument('--ylim', type=float, nargs=2)



        args = parser.parse_args()


    else:
        os.chdir(os.getcwd().replace("manuscript_figure_script", ""))
        parser = argparse.ArgumentParser()
        args = parser.parse_args()
        args.simul_data  = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_segsites.tsv"
        args.llk_rand    = "simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_gridsearch_R0"

        args.figname = "figure2_R0=16e-1_seed1234"
        args.ylim = [-90, -72.5]
        args.xlim = [1, 2.5]
        args.true_param  = [1.6]


    print(args)
    hey = run(args.simul_data, args.llk_rand, args.true_param, args.xlim, args.ylim,  args.figname)


    print(type(args.xlim[1]))










