import os
import numpy as np
import argparse

import pandas as pd
from matplotlib import pyplot as plt

def check_same_npz_easy(dir_npz1, dir_npz2):

    npz1 = np.load(dir_npz1, allow_pickle=True)
    npz2 = np.load(dir_npz2, allow_pickle=True)

    list_npz = list(set(npz1.files).intersection(set(npz2.files)))

    print(list_npz)

    diff = {}
    for file in list_npz:
        if not np.array_equal(npz1[file], npz2[file], equal_nan = True):
            #print(npz1[file], npz2[file])
            print(f"{file} is different!")
            diff[file] =[npz1[file], npz2[file]]

    if len(diff) > 0:
        print([x[0] for x in diff])
        for i in diff:
            print(i[0])
            print(i[1], i[2])
            print("- - - -")
        #print(npz1[file], npz2[file])
        #raise Exception(f"{diff[0]} is different!")
    else:
        print(">> two npzs are same!")

    return






def check_same_npz(dir_npz1, dir_npz2):

    npz1 = np.load(dir_npz1, allow_pickle=True)
    npz2 = np.load(dir_npz2, allow_pickle=True)

    if npz1.files != npz2.files:
        print(npz1.files)
        print(npz2.files)
        raise Exception("list of files are different!")

    diff = None
    for i, file in enumerate(npz1.files):
        print(f'{i}/{len(npz1.files)}: {file}')
        if not np.array_equal(npz1[file], npz2[file], equal_nan = True):
            print(f"{file} is different!")
            diff =  file, npz1[file], npz2[file]

    if not diff is None:
        raise Exception(f"{diff[0]} is different!")
    else:
        print(">> two npzs are same!")

    return

def check_same_statevar(dir_tsv1, dir_tsv2):
    a = pd.read_csv(dir_tsv1).to_numpy()
    b = pd.read_csv(dir_tsv2).to_numpy()

    if np.array_equal(a, b):
        print(">> two tsv are same!")
    else:
        raise Exception(">> two tsv are different!")

    return

def find_other_trials(fname):
    import re

    var_part = r'_trial[0-9]+'
    match = re.search(var_part, fname)

    if match:
        parent_folder = fname[:match.span()[0]].rpartition("/")[0]
        var_names = fname[:match.span()[0]] + var_part + fname[match.span()[1]:]
        var_names = var_names.replace(parent_folder, "")

        print(var_names)
        return [x[0]+"/"+y for x in os.walk(parent_folder) for y in x[2] if re.search(var_names, x[0]+"/"+y)]
        #return [f'{fname_folder}/{x}' for x in os.listdir(fname_folder) if re.search(var_names, x)]
    else:
        return [fname]

def check_simpleSEIR_hSEIR(dir1, dir2):
    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(4, 1, figsize=[4, 6], dpi=300)

    def plot_statevar_hSEIR(ax, label, c, filename, time_shift=0, ls="-", alpha=0.7, compartment = "I"):
        statevar = pd.read_csv(filename, sep="\t", comment="#", names=["t", "S", "E", "I_l", "I_h", "R"])
        if compartment == "I":
            ax.plot(statevar['t'] + time_shift, statevar['I_l'] + statevar['I_h'], alpha=alpha, label=label, c=c, ls=ls)
        else:
            ax.plot(statevar['t'] + time_shift, statevar[compartment], alpha=alpha, label=label, c=c, ls=ls)
        ax.set_xlabel('Day post index case')
        ax.set_ylabel(f'# of {compartment}s')

    def plot_statevar_simpleSEIR(ax, label, c, filename, time_shift=0, ls="-", alpha=0.7, compartment = "I"):
        statevar = pd.read_csv(filename, sep="\t", comment="#", names=["t", "S", "E", "I","R"])
        ax.plot(statevar['t'] + time_shift, statevar[compartment], alpha=alpha, label=label, c=c, ls=ls)

        ax.set_xlabel('Day post index case')
        ax.set_ylabel(f'# of {compartment}s')


    label = "simpleSEIR"
    for file in find_other_trials(dir1):
        plot_statevar_simpleSEIR(axes[0], label, "r", file, alpha=0.3, compartment="S")
        plot_statevar_simpleSEIR(axes[1], label, "r", file, alpha=0.3, compartment="E")
        plot_statevar_simpleSEIR(axes[2], label, "r", file, alpha=0.3, compartment="I")
        plot_statevar_simpleSEIR(axes[3], label, "r", file, alpha=0.3, compartment="R")
        label = None

    label = "SEI1I2R"
    for file in find_other_trials(dir2):
        plot_statevar_hSEIR(axes[0], label, "b", file, alpha = 0.3, compartment="S")
        plot_statevar_hSEIR(axes[1], label, "b", file, alpha=0.3, compartment="E")
        plot_statevar_hSEIR(axes[2], label, "b", file, alpha=0.3, compartment="I")
        plot_statevar_hSEIR(axes[3], label, "b", file, alpha=0.3, compartment="R")
        label = None


    axes[0].axvline(x=24, linewidth=0.5, ls='--')
    axes[0].axvline(x=28, linewidth=0.5, ls='--')
    axes[0].axvline(x=32, linewidth=0.5, ls='--')
    axes[0].legend()
    fig.tight_layout(pad = 0.5)

    os.makedirs('simulated_data/test', exist_ok=True)
    fig.savefig('simulated_data/test/check_simpleSEIR_hSEIR.png', dpi =300)
    return

def n_segregating_over_time(ax, data, c = "k", label = None, alpha=1):
    data = pd.read_csv(data, sep = "\t")
    data = data.loc[~np.isnan(data['s'])]
    #print(data)
    #print("\n\n\n\n")

    ax.plot   (data['window_end'], data['s'], c=c, label=label, alpha=alpha)
    ax.scatter(data['window_end'], data['s'], c=c, label=label, alpha=alpha, s= 15)

    ax.set_ylim([0, int(np.nanmax(data['s']) * 1.1)])
    ax.set_xlabel("Days post index case")
    ax.set_ylabel("Number of segregating sites")
    return

def n_sequence_over_time(ax, data, c = "k", label = None, alpha=1):
    data = pd.read_csv(data, sep = "\t")
    data = data.loc[~np.isnan(data['n'])]
    ax.plot   (data['window_end'], data['n'], c=c, label=label, alpha=alpha)
    ax.scatter(data['window_end'], data['n'], c=c, label=label, alpha=alpha, s=15)

    #ax.set_xlabel("Days post index case")
    ax.set_ylabel("Number of sequences")

    return

def sampling_prop_over_time(ax, data, c = "k", label = None, alpha=1):
    data = pd.read_csv(data, sep = "\t")
    data = data.loc[~np.isnan(data['n'])]
    #print(data)
    #print("\n\n\n\n")

    ax.plot   (data['window_end'], data['n'] / data['n_recovered_in_win'], c=c, label=label, alpha=alpha, linestyle = "--", linewidth = 0.7)
    ax.scatter(data['window_end'], data['n'] / data['n_recovered_in_win'], c=c, label=label, alpha=alpha, s=8, marker = "s")

    #ax.set_xlabel("Days post index case")
    ax.set_ylabel("Sampling proportion")

    return



def n_recovered_over_time(ax, data, c="k", label=None, alpha=1):
        data = pd.read_csv(data, sep="\t")
        data = data.loc[~np.isnan(data['n'])]
        ax.plot(data['window_end'], data['n_recovered_in_win'], c=c, label=label, alpha=alpha)
        ax.scatter(data['window_end'], data['n_recovered_in_win'], c=c, label=label, alpha=alpha, s=15)

        # ax.set_xlabel("Days post index case")
        ax.set_ylabel("Number of recovered")

        return


def check_simpleSEIR_hSEIR_segsites(dir1, dir2):
    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(1, 1, figsize=[3, 4], dpi=300)


    label = "simpleSEIR"
    for file in find_other_trials(dir1):
        n_segregating_over_time(axes, file, label=label, c="r", alpha=0.3)
        label = None

    label = "SEI1I2R"
    for file in find_other_trials(dir2):
        n_segregating_over_time(axes, file, label=label, c="b", alpha=0.3)
        label = None


    axes.legend()
    fig.tight_layout(pad = 0.5)

    os.makedirs('simulated_data/test', exist_ok=True)
    fig.savefig('simulated_data/test/check_simpleSEIR_hSEIR_segsites.png', dpi =300)
    return


def show_segsite_trajctories(dir1):
    plt.rcParams.update({'font.size': 7})
    fig, axes = plt.subplots(4, 1, figsize=[4, 8], dpi=300)

    label = None
    for file in find_other_trials(dir1):
        n_recovered_over_time(axes[0], data=dir1)
        sampling_prop_over_time(axes[1], data=dir1)
        n_sequence_over_time    (axes[2], data=dir1)
        n_segregating_over_time (axes[3], data=dir1)
        label = None


    fig.suptitle("\n".join(dir1.split("/")))
    fig.tight_layout()
    os.makedirs('simulated_data/test', exist_ok=True)
    fig.savefig(f'{dir1.rpartition("/")[0]}/check_segsite_trajectory.png', dpi =300)
    return


def show_epidynamics(dir1):
    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(4, 1, figsize=[4, 6], dpi=300)

    def plot_statevar_simpleSEIR(ax, label, c, filename, time_shift=0, ls="-", alpha=0.7, compartment = "I"):
        statevar = pd.read_csv(filename, sep="\t", comment="#", names=["t", "S", "E", "I","R"])
        ax.plot(statevar['t'] + time_shift, statevar[compartment], alpha=alpha, label=label, c=c, ls=ls)

        ax.set_xlabel('Day post index case')
        ax.set_ylabel(f'# of {compartment}s')

    label = None
    for file in find_other_trials(dir1):
        plot_statevar_simpleSEIR(axes[0], label, "r", file, alpha=0.3, compartment="S")
        plot_statevar_simpleSEIR(axes[1], label, "r", file, alpha=0.3, compartment="E")
        plot_statevar_simpleSEIR(axes[2], label, "r", file, alpha=0.3, compartment="I")
        plot_statevar_simpleSEIR(axes[3], label, "r", file, alpha=0.3, compartment="R")
        label = None

    fig.tight_layout(pad = 0.5)

    os.makedirs('simulated_data/test', exist_ok=True)
    fig.savefig(f'{dir1.rpartition("/")[0]}/{dir1.rpartition("/")[-1].replace(".tsv", ".png")}', dpi =300)
    #fig.savefig('simulated_data/test/check_simpleSEIR_hSEIR.png', dpi =300)
    return



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(dest='test')

    args_check_same_npz = subparser.add_parser('check_same_npz')
    args_check_same_npz.add_argument('dir1')
    args_check_same_npz.add_argument('dir2')


    args_check_same_npz_easy = subparser.add_parser('check_same_npz_easy')
    args_check_same_npz_easy.add_argument('dir1')
    args_check_same_npz_easy.add_argument('dir2')


    args_check_same_statevar_npz = subparser.add_parser('check_same_statevar')
    args_check_same_statevar_npz.add_argument('dir1')
    args_check_same_statevar_npz.add_argument('dir2')

    args_check_simpleSEIR_hSEIR = subparser.add_parser('check_simpleSEIR_hSEIR')
    args_check_simpleSEIR_hSEIR.add_argument('dir1')
    args_check_simpleSEIR_hSEIR.add_argument('dir2')

    args_check_simpleSEIR_hSEIR_segsites = subparser.add_parser('check_simpleSEIR_hSEIR_segsites')
    args_check_simpleSEIR_hSEIR_segsites.add_argument('dir1')
    args_check_simpleSEIR_hSEIR_segsites.add_argument('dir2')


    args_check_segsite_trajctories = subparser.add_parser('check_segsite_trajctories')
    args_check_segsite_trajctories.add_argument('dir1')

    args_show_epidynamics = subparser.add_parser('show_epidynamics')
    args_show_epidynamics.add_argument('dir1')



    args = parser.parse_args()

    if args.test == 'check_same_npz':
        check_same_npz(args.dir1, args.dir2)
    if args.test == 'check_same_npz_easy':
        check_same_npz_easy(args.dir1, args.dir2)
    elif args.test == 'check_same_statevar':
        check_same_statevar(args.dir1, args.dir2)
    elif args.test == 'check_simpleSEIR_hSEIR':
        check_simpleSEIR_hSEIR(args.dir1, args.dir2)
    elif args.test == 'check_simpleSEIR_hSEIR_segsites':
        check_simpleSEIR_hSEIR_segsites(args.dir1, args.dir2)
    elif args.test == 'check_segsite_trajctories':
        show_segsite_trajctories(args.dir1)
    elif args.test == 'show_epidynamics':
        show_epidynamics(args.dir1)


#
#     args.test = "check_same_npz_easy"
#     args.dir1 = "/Users/yeongseon/Dropbox/PhD_Emory_university/1_manuscripts/(manuscript) segregating_site/jey/manuscript_figures_data/fig4/SEIR_simple_simtaj.npz"
#     args.dir2 = "/Users/yeongseon/Dropbox/PhD_Emory_university/1_projects/Project_seg_site/simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_simtaj.npz"
#
#     check_same_npz_easy(args.dir1, args.dir2)
#
# ## genotype_info
# trim1 = np.where(diff['genotype_info'][0].sum(axis=1) > 0)[0].max()
# trim2 = np.where(diff['genotype_info'][1].sum(axis=1) > 0)[0].max()
# np.array_equal(diff['genotype_info'][0][:trim1], diff['genotype_info'][1][:trim2], equal_nan=True)
#
# ## n_segregating contains no meaning in this simulation
#
# ## recovered
# np.array_equal(diff['recovered'][0][:, :2], diff['recovered'][1][:, :2], equal_nan=True)
#
# diff['statevar'][0][:, -1] -= 1
# np.array_equal(diff['statevar'][0], diff['statevar'][1], equal_nan=True)
