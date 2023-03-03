import os
import pickle
import datetime
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator

from fig4_plot import n_segregating_over_time, plot_state_variable
from fig7_plot import convert_to_matlab_date, convert_to_caldate
# ---------------------------
def format_func(x, tick_number):
    return convert_to_caldate(x).replace("2020/\n", "").replace("2019/\n", "")

def extract_statevar_and_segsite(dir, n_patricle_to_extract):
    reconstructed_files = [x for x in os.listdir(dir) if (x.find(").pkl") > -1)]

    statevar = []
    infected_outside = []
    n_segregating = []

    for reconstructed_file in reconstructed_files:
        particles = pickle.load(open(dir + "/" + reconstructed_file, "rb"))
        random_particle_idx = np.random.choice(len(particles), size=n_patricle_to_extract, replace=False)
        # print ("reconstructed_file:", reconstructed_file)
        for particle_idx in random_particle_idx:
            particle = particles[particle_idx]

            trim_statevar = particle[0][:, 0:-2].sum(axis=1) > 0
            trim_statevar = (convert_to_matlab_date("2019/12/24") <= particle[0][:, 0])
            statevar.append(particle[0][trim_statevar])

            print(statevar[-1][:5, :])

            tmp_infected_outside = np.array(particle[1])
            trim_infected = (convert_to_matlab_date("2019/12/24") <= tmp_infected_outside[:, 0])
            infected_outside.append(tmp_infected_outside[trim_infected])

            trim_s = (737839 <= particle[1][:, 1]) * (particle[1][:, 1] <= 737867)
            n_segregating.append(particle[1][trim_s])

    return statevar, n_segregating, infected_outside





def run(args):
    txt = open(f'manuscript figures/txt/{args.figname}.txt', 'w')
    plt.rcParams.update({'font.size': 10})
    fig, axes =  plt.subplots(4, 3, figsize=[8, 9], dpi=300)
    # ---------------------------
    dir_observed_dataset = args.observed_data

    ## reconstructed state variables
    for i, dir_particle in enumerate(args.dir_particle):
        df_particle_statevar, df_n_segregating, infected_outside = extract_statevar_and_segsite(dir_particle, 1)  ## n_patricle_to_extract_per_combo

        compartments = ['S', 'E1', 'E2_l', 'I_l', 'E2_h', 'I_h', 'R', 'O_l', 'O_h']
        compartments = {compartments[i]: i - 1 for i in range(len(compartments))}
        idx_all_infected = np.array([compartments[x] for x in ['E1', 'E2_l', 'I_l', 'E2_h', 'I_h']])
        idx_sum_for_N = np.array([compartments[x] for x in ['S', 'E1', 'E2_l', 'I_l', 'E2_h', 'I_h', 'R']])
        idx_sum_for_cumR = np.array([compartments[x] for x in ['E1', 'E2_l', 'I_l', 'E2_h', 'I_h',  'R']])


        cum_R = []
        total_infected = []
        for j in np.random.choice(len(df_n_segregating), size=10):
            n_segregating_over_time(axes[0, i], df_n_segregating[j], particle=True)

            df_state = df_particle_statevar[j]

            denom       =         df_state[0, idx_sum_for_N   +2].sum()
            state_cumR  =        ((df_state[:,idx_sum_for_cumR+2]/denom) * 100).sum(axis=1)
            state_all_infected = (df_state[:,idx_all_infected +2].sum(axis=1)/denom) * 100

            axes[1, i].plot(df_state[:, 0], state_all_infected, c="k", alpha=0.1, label="_nolegend_")
            axes[2, i].plot(df_state[:, 0], state_cumR, c="k", alpha=0.1, label="_nolegend_")
            axes[3, i].plot(infected_outside[j][:, 0], infected_outside[j][:, 1].cumsum(), c="k", alpha=0.1, label="_nolegend_")

            # print (convert_to_caldate(df_state[:, 0][state_all_infected>0][0]).replace("2020/\n", "").replace("2019/\n", ""),
            #        df_state[:, 2:-1].sum(axis=1)[state_all_infected>0][0],
            #        convert_to_caldate(infected_outside[j][:, 0][infected_outside[j][:,1].cumsum() > 0][0]).replace("2020/\n", "").replace("2019/\n", ""),
            #        infected_outside[j][:, 1].cumsum()[infected_outside[j][:, 1].cumsum() > 0][0])

            ## save for printing
            filter = (df_state[:, 0] < 737864 + 0.0001) * ( 737864 - 0.0001 < df_state[:, 0])
            total_infected.append(state_all_infected[filter][0])
            cum_R.append(state_cumR[filter][0])



        ## Observed segregating site & serology data
        df_simulated_dataset = pd.read_csv(dir_observed_dataset, sep="\t")

        n_segregating_over_time(axes[0, i], df_simulated_dataset, particle=False)
        axes[2, i].errorbar([737864], [0.0041*100], yerr=[[0.0041*100 - 0.0005*100], [0.0088*100 - 0.0041*100]], color="black")
        axes[2, i].scatter ([737864], [0.0041*100], color="black", s=7)

        ## print
        print(f"{min(cum_R):.3f}%, {max(cum_R):.3f}", file=txt)
        print(f"{min(cum_R):.3f}%, {max(cum_R):.3f}")
        print ("-----", file=txt)

        ## set axis ticks & axis tick labels
        xticks_segregating = df_simulated_dataset.loc[~np.isnan(df_simulated_dataset['s'])]['window_end']
        x_ticks_statevar = np.arange(xticks_segregating.max(), convert_to_matlab_date("2019/12/23"), -14)[::-1]
        for j in range(4):
            if j == 0:
                #axes[0, i].xaxis.set_major_locator(FixedLocator(xticks_segregating))  # FixedLocator
                axes[0, i].set_xlim([x_ticks_statevar.min(), x_ticks_statevar.max()+2])
                axes[0, i].xaxis.set_major_locator(FixedLocator(x_ticks_statevar))  # FixedLocator

            else:
                axes[j, i].set_xlim([x_ticks_statevar.min(), x_ticks_statevar.max()+2])
                axes[j, i].xaxis.set_major_locator(FixedLocator(x_ticks_statevar))  # FixedLocator
                if j != 3:
                    axes[j, i].set_yscale('log')

            if j == 3:
                axes[j, i].yaxis.set_major_locator(FixedLocator([1e0, 1e1, 1e2, 1e3, 1e4, 1e5]))

            axes[j, i].xaxis.set_major_formatter(plt.FuncFormatter(format_func))
            plt.setp(axes[j, i].get_xticklabels(), rotation=30, ha='right', rotation_mode='anchor')

        axes[0, 0].set_ylabel("Number of\nsegregating sites", fontsize=9)
        axes[1, 0].set_ylabel("Infected individuals (%)", fontsize=9)
        axes[2, 0].set_ylabel("Seroprevalence (%)", fontsize=9)
        axes[3, 0].set_ylabel("Cumulative infections from\noutside-of-France contacts", fontsize=8)

        axes[0, i].set_xlabel(None, fontsize=9)
        axes[3, i].set_xlabel("Days post index case", fontsize=9)

        axes[1, i].set_ylim([axes[1, i].get_ylim()[0], 0.015 * 100])
        axes[2, i].set_ylim([axes[2, i].get_ylim()[0], 0.015 * 100])

    # ---------------------------
    ## add subplot label
    labels = ["ABC".lower(), "DEF".lower(), "GHI".lower(), "JKL".lower()]
    for i, label_row in enumerate(labels):
        for j, l in enumerate(label_row):
            if j == 0:
                x_loc = -0.4
                axes[i][j].text(-0.4, 1.1, l, transform=axes[i][j].transAxes, fontsize=16, fontweight='bold', va='top',
                                ha='left')
            else:
                x_loc = -0.3
                axes[i][j].text(-0.3, 1.1, l, transform=axes[i][j].transAxes, fontsize=16, fontweight='bold', va='top',
                                ha='left')
            continue

    ## save figure
    fig.tight_layout()
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
        parser.add_argument('--figname', type=str, default="figure9")
        parser.add_argument('--observed_data', type=str, required=True)
        parser.add_argument('--dir_particle' , type=str, nargs = 3, required=True)

        args = parser.parse_args()


    else:
        os.chdir(os.getcwd().replace("manuscript_figure_script", ""))
        parser = argparse.ArgumentParser()
        args = parser.parse_args()

    print(args)
    hey = run(args)

