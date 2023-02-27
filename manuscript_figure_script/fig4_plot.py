import os
import pickle
import datetime
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator
# ---------------------------
def n_segregating_over_time(axes, data, particle=False):
    if particle:
        axes.plot(data[:, 1], data[:, 0], c="k", alpha=0.15, linewidth=1, label="_nolegend_")
    else:
        data = data.loc[~np.isnan(data['s'])]
        axes.plot(data['window_end'], data['s'], linewidth=1.25, c="r", linestyle='dashed')
        # axes.scatter(data['window_end'], data['s'], linewidth=1.5, c="r", s=6)

    axes.set_xlabel("Days post index case", fontsize=9)
    axes.set_ylabel("Number of\nsegregating sites", fontsize=9)

    axes.xaxis.set_label_coords(.45, -.2)

    return


def plot_state_variable(axes, df_state, statavar_no, ylabel, particle=False, ylim=None):
    if particle:
        axes.plot(df_state[:, 0], df_state[:, statavar_no], c="k", alpha=0.15, linewidth=1, label="_nolegend_")
    else:
        axes.plot(df_state[:, 0], df_state[:, statavar_no], c="r", linewidth=1.25, linestyle='dashed')

    axes.set_xlabel("Days post index case", fontsize=9)
    axes.set_ylabel("Number of\n" + ylabel + " individuals", fontsize=9)

    axes.xaxis.set_label_coords(.45, -.2)

    if ylim:
        axes.set_ylim(ylim)

    return

def extract_statevar_and_segsite(dir, n_patricle_to_extract):
    reconstructed_files = [x for x in os.listdir(dir) if (x.find(".pkl") > -1) and (x.find("(") > -1)]

    statevar = []
    n_segregating = []

    for reconstructed_file in reconstructed_files:
        particles = pickle.load(open(dir + "/" + reconstructed_file, "rb"))
        random_particle_idx = np.random.choice(len(particles), size=n_patricle_to_extract, replace=False)
        # print ("reconstructed_file:", reconstructed_file)
        for particle_idx in random_particle_idx:
            particle = particles[particle_idx]

            trim_statevar = particle[0][:, 1:4].sum(axis=1) > 0
            statevar.append(particle[0][trim_statevar])

            trim_s = (0 <= particle[1][:, 1]) * (particle[1][:, 1] <= 154.0)
            n_segregating.append(particle[1][trim_s])

    return statevar, n_segregating



def run(args):

    dir_simulated_dataset = args.true_segregating
    dir_true_statevar     = args.true_sim
    dir_particle          = args.reconstructed

    df_simulated_dataset = pd.read_csv(dir_simulated_dataset, sep = "\t")
    df_true_statevar     = np.load    (dir_true_statevar, allow_pickle = True)["statevar"]
    # ---------------------------
    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(1, 4, figsize=[8, 2.25], dpi=300)

    ## (A, B, C, D) reconstructed state variables
    particle = True
    df_particle_statevar, df_n_segregating = extract_statevar_and_segsite(dir_particle, 1)           ## n_patricle_to_extract_per_combo
    for i in np.random.choice(len(df_n_segregating), size = 30):
        n_segregating_over_time(axes[0], df_n_segregating[i], particle=True)
        plot_state_variable(axes[1], df_particle_statevar[i], 1, 'susceptible', particle=True)   # [0, 110000])
        plot_state_variable(axes[2], df_particle_statevar[i], 2, 'exposed', particle=True)       # [0, 6000])
        plot_state_variable(axes[3], df_particle_statevar[i], 3, 'infected', particle=True)      # [0, 8000])

    n_segregating_over_time(axes[0], df_simulated_dataset, particle=False)
    plot_state_variable    (axes[1], df_true_statevar, 1, 'susceptible', particle=False)    #  [0, 110000])
    plot_state_variable    (axes[2], df_true_statevar, 2, 'exposed',  particle=False)        #  [0, 6000])
    plot_state_variable    (axes[3], df_true_statevar, 3, 'infected',  particle=False)       #  [0, 8000])

    # ---------------------------
    # axis
    axes[1].ticklabel_format(axis='y', style="sci", scilimits=(3, 4))
    axes[1].yaxis.major.formatter._useMathText = True

    ## add subplot label
    labels = 'ABCD'.lower()
    for i, label_row in enumerate(labels):
        axes[i].text(-0.55, 1.15, label_row, transform=axes[i].transAxes, fontsize=16, fontweight='bold', va='top',
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
        parser.add_argument('--figname', type=str, default="figure3")
        parser.add_argument('--true_sim', type=str, required=True)
        parser.add_argument('--true_segregating', type=str, required=True)
        parser.add_argument('--reconstructed', type=str, required=True)



        args = parser.parse_args()


    else:
        os.chdir(os.getcwd().replace("manuscript_figure_script", ""))
        parser = argparse.ArgumentParser()
        args = parser.parse_args()



    print(args)
    hey = run(args)


