import argparse
import baltic as bt
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
from utils import datetime_from_numeric
import string
import pandas as pd



def plot_divergence_tree(ax, tree_file=None, tree_name_sep="|", tree_name_field=1, metadata_dict={}, plot_config={}):    
    myTree = bt.loadNewick(tree_file)
    myTree.sortBranches()
    myTree.drawTree()
    xMin = min([i.x for i in myTree.Objects])
    xMax = max([i.x for i in myTree.Objects])
    xSpan = xMax - xMin
    myTree.plotTree(ax, 
        x_attr=lambda k: k.x,
        width=plot_config['branch_width'] if 'branch_width' in plot_config.keys() else 1.0, 
        colour=plot_config['colors']['base'])
    myTree.plotPoints(ax, 
        x_attr=lambda k: k.x, 
        size=10, 
        colour=lambda k: plot_config['colors'][metadata_dict[k.name.split(tree_name_sep)[tree_name_field]]],
        outline_colour=plot_config['colors']['base'])
    ax.set_xlim((xMin-xSpan*0.1, xMax+xSpan*0.1))
    ax.set_ylim(-myTree.ySpan*0.01, myTree.ySpan*1.01)
    ax.set_yticks([])
    x_tick_labels = np.arange(np.ceil(xMin*29903-xMin*29903), np.ceil(xMax*29903-xMin*29903)).astype(int)
    x_ticks = x_tick_labels/29903 + xMin
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)
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




def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--trees', 
        default=None, help='newick format trees to plot, takes two', nargs=2)
    parser.add_argument('--treeNameSep', default='|')
    parser.add_argument('--metadata')
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

    fig, axs = plt.subplots(1,2, figsize=(6.4*1.25, 4.8*1.5), constrained_layout=True)
    axs[0] = plot_divergence_tree(axs[0], tree_file=args.trees[0], tree_name_sep=args.treeNameSep, tree_name_field=args.treeNameField, metadata_dict=metadata_dict, plot_config=plot_config)
    axs[1] = plot_time_tree(axs[1], tree_file=args.trees[1], tree_name_sep=args.treeNameSep, tree_name_field=args.treeNameField, metadata_dict=metadata_dict, plot_config=plot_config)
    for ax_idx, ax in enumerate(axs):
            ax.text(0, 1.05, string.ascii_uppercase[ax_idx], transform=ax.transAxes, 
                size=16, weight='bold', va="top", color='#333333')

    fig.savefig(args.outName + '.pdf')
    plt.close()



if __name__ == "__main__":
    run()

