import os
import argparse
import baltic as bt
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
from utils import datetime_from_numeric
import string
import pandas as pd
from scipy.stats import linregress
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib.dates import MonthLocator, DateFormatter
import matplotlib.dates as mdates


def numeric_from_datetime(dt):
    from calendar import isleap
    import datetime
    if dt is None:
        dt = datetime.datetime.now()
    days = 366 if isleap(dt.year) else 365
    # removing the 0.5 adjustment in the standard treetime code 
    # to align with tajimas D inference method
    #res = dt.year + (dt.timetuple().tm_yday-0.5) / days
    res =  dt.year + (dt.timetuple().tm_yday) / days
    return(res)




def plot_divergence_tree(ax, tree_file=None, tree_name_sep="|",
 tree_name_field=1, metadata_dict={}, plot_config={},
 tree_length_multiplier=1):    
    myTree = bt.loadNewick(tree_file)
    myTree.sortBranches()
    myTree.drawTree()
    xMin = min([i.x for i in myTree.Objects])
    xMax = max([i.x for i in myTree.Objects])
    xSpan = (xMax - xMin) * tree_length_multiplier
    myTree.plotTree(ax, 
        x_attr=lambda k: (k.x - xMin)*tree_length_multiplier if k.x else 0,
        width=plot_config['branch_width'] if 'branch_width' in plot_config.keys() else 1.0, 
        colour=plot_config['colors']['base'])
    myTree.plotPoints(ax, 
        x_attr=lambda k: (k.x - xMin)*tree_length_multiplier if k.x else 0, 
        size=10, 
        colour=lambda k: plot_config['colors'][metadata_dict[k.name.split(tree_name_sep)[tree_name_field]]],
        outline_colour=plot_config['colors']['base'])
    ax.set_xlim(-xSpan*0.1, xSpan*1.1)
    ax.set_ylim(-myTree.ySpan*0.01, myTree.ySpan*1.01)
    ax.set_yticks([])
    #x_tick_labels = np.arange(np.ceil(xMin*29903-xMin*29903), np.ceil(xMax*29903-xMin*29903)).astype(int)
    #x_ticks = x_tick_labels/29903 + xMin
    #ax.set_xticks(x_ticks)
    #ax.set_xticklabels(x_tick_labels)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('substitutions')
    [ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'left']]
    ax.text(0.5, 0.65, 'France (included)', transform=ax.transAxes, 
            size=16, va="top", color=plot_config['colors']['included'])
    ax.text(0.5, 0.61, 'France (excluded)', transform=ax.transAxes, 
            size=16, va="top", color=plot_config['colors']['not_included_france'])
    ax.text(0.5, 0.57, 'Other (excluded)', transform=ax.transAxes, 
            size=16, va="top", color=plot_config['colors']['not_included_exog'])
    return(ax)



def plot_time_tree(ax, tree_file=None, tree_name_sep="|", tree_name_field=1, metadata_dict={}, plot_config=None):
    myTree = bt.loadNewick(tree_file)
    myTree.sortBranches()
    myTree.drawTree()
    myTree.setAbsoluteTime(plot_config['max_date'])
    xMin = min([i.absoluteTime for i in myTree.Objects])
    xMax = max([i.absoluteTime for i in myTree.Objects])
    xSpan = xMax - xMin
    L=len(list(filter(lambda k: k.branchType=='leaf',myTree.Objects)))
    tip_size = 50+myTree.treeHeight*35
    myTree.plotTree(ax, 
        x_attr=lambda k: k.absoluteTime,
        width=1.0, 
        colour=plot_config['colors']['base'])
    ax.set_xlim((xMin-xSpan*0.1, xMax+xSpan*0.1))
    myTree.plotPoints(ax, 
        x_attr=lambda k: k.absoluteTime, 
        size=10, 
        colour=lambda k: plot_config['colors'][metadata_dict[k.name.split(tree_name_sep)[tree_name_field]]],
        outline_colour=plot_config['colors']['base'])
    ax.set_ylim(-myTree.ySpan*0.01, myTree.ySpan*1.01)
    ax.set_yticks([])
    [ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'left']]
    x_ticks = ax.get_xticks()
    # todo make this not manual and lazy
    x_ticks = [2020+1/366, 2020+15/366, 2020+(31+1)/366, 2020+(31+15)/366, 2020+(31+29+1)/366, 2020+(31+29+15)/366]
    x_tick_labels = [datetime_from_numeric(i).strftime("%Y-%m-%d") for i in x_ticks]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)
    _ = [tick.set_rotation(45) for tick in ax.get_xticklabels()]
    _ = [tick.set_horizontalalignment('right') for tick in ax.get_xticklabels()]
    ax.set_xlabel('date')
    return(ax)




def get_dist_dat(dist_path, date_dict):
    dist_dat = pd.read_csv(dist_path, header=None, sep='\t')
    dist_dat[dist_dat.shape[1]] = dist_dat[0].apply(lambda k: date_dict[k.split('|')[1]])
    # regress on the dist dat
    regress = linregress(dist_dat[dist_dat.shape[1]-1].apply(numeric_from_datetime), 
        dist_dat[1])
    return(dist_dat, regress)



def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--distDat', default=None)
    parser.add_argument('--trees', 
        default=None, help='newick format trees to plot, takes two', nargs=2)
    parser.add_argument('--treeNameSep', default='|')
    parser.add_argument('--metadata')
    parser.add_argument('--treeLengthMultiplier', type=int, default=29903)
    parser.add_argument('--treeNameField',
        default=1, type=int, 
        help='which field in the sequence name corresponds to the metadata name column')
    parser.add_argument('--outName', default='trees')
    args = parser.parse_args()

    #args.metadata = 'data/france_random_global_ref_figure.tsv' 
    #args.trees = \
    #    ['data/france_random_global_ref_aligned_refine.nwk',
    #        'data/france_random_global_ref_aligned_refine_time.nwk']
    #args.outName = 'figures/france_trees'

    plot_config = {"colors": 
        {"base": "#575c66",
            "included": "#0D7BCA",
            "not_included_france": "#9ecae9",
            "not_included_exog": "#575c66"}, 
        'max_date': 2020 + (31+29+17)/366}


    metadata = pd.read_csv(args.metadata, sep='\t', header=None)
    metadata_dict = {i[0]:i[1] for i in metadata.values}

    metadata[2] = pd.to_datetime(metadata[2])
    date_dict = {i[0]: i[2] for i in metadata.values}

    dist_dat, regress = \
        get_dist_dat(args.distDat, date_dict)

    #fig, axs = plt.subplots(1,2, figsize=(6.4*1.25, 4.8*1.5), constrained_layout=True)
    fig = plt.figure(figsize=(6.4*1.2, 4.8*2), constrained_layout=True)
    gs = gridspec.GridSpec(6, 2, figure=fig)

    axs = [plt.subplot(gs[:4,0]), 
        plt.subplot(gs[:4,1]), 
        plt.subplot(gs[4:,:])]


    axs[0] = plot_divergence_tree(axs[0], tree_file=args.trees[0], 
        tree_name_sep=args.treeNameSep, tree_name_field=args.treeNameField, 
        metadata_dict=metadata_dict, plot_config=plot_config,
        tree_length_multiplier=args.treeLengthMultiplier)
    axs[1] = plot_time_tree(axs[1], tree_file=args.trees[1], tree_name_sep=args.treeNameSep, tree_name_field=args.treeNameField, metadata_dict=metadata_dict, plot_config=plot_config)
    regress_range = \
            (dist_dat[dist_dat.shape[1]-1].min(), 
                dist_dat[dist_dat.shape[1]-1].max())
    axs[2].scatter(dist_dat[2], dist_dat[1], 
        alpha=0.25, color='#81a1c1')
    axs[2].plot(regress_range, 
        [regress.slope*(numeric_from_datetime(x)) + 
            regress.intercept for x in regress_range], 
        color='#BF616A')
    axs[2].text(0.05, 0.8, 
        str(round(regress.slope)) + " subs/yr, \n{:.2e} subs/site/yr".format(regress.slope/29903),
        transform=axs[2].transAxes, color='#BF616A', size=12)
    axs[2].set_ylabel('distance to Wuhan/Hu-1')
    axs[2].xaxis.set_major_locator(mdates.DayLocator(interval=5))
    #days = mdates.DayLocator(interval=7)
    #axs[2].xaxis.set_minor_locator(days)
    #axs[2].xaxis.set_minor_formatter(DateFormatter("%Y-%m-%d"))
    axs[2].xaxis.set_major_formatter(DateFormatter("%Y-%m-%d"))
    _ = [tick.set_rotation(45) for tick in axs[2].get_xticklabels()]
    _ = [tick.set_horizontalalignment('right') for tick in axs[2].get_xticklabels()]
    x_pos = [0,0,-0.075]
    for ax_idx, ax in enumerate(axs):
            ax.text(x_pos[ax_idx], 1.05, string.ascii_lowercase[ax_idx], transform=ax.transAxes,
                size=16, weight='bold', va="top", color='#333333')

    try:
        fig.savefig(f'{args.outName.rpartition("/")[0]}/pdf/{args.outName.rpartition("/")[-1]}.pdf')
        fig.savefig(f'{args.outName.rpartition("/")[0]}/png/{args.outName.rpartition("/")[-1]}.png')
    except:
        os.makedirs(f'{args.outName.rpartition("/")[0]}/pdf/{args.outName.rpartition("/")[-1]}.pdf')
        os.makedirs(f'{args.outName.rpartition("/")[0]}/png/{args.outName.rpartition("/")[-1]}.png')
        fig.savefig(f'{args.outName.rpartition("/")[0]}/pdf/{args.outName.rpartition("/")[-1]}.pdf')
        fig.savefig(f'{args.outName.rpartition("/")[0]}/png/{args.outName.rpartition("/")[-1]}.png')

    plt.close()



if __name__ == "__main__":
    run()

