import argparse
import subprocess
import shlex
from Bio import SeqIO, Phylo
import Bio
from treetime import TreeTime
from datetime import date
import pandas as pd
import numpy as np


'''
def infer_phylogeny(aln_path, n_bootstrap, prefix):
    if n_bootstrap != 'None':
        iqtree_cmd = iqtree_cmd = f'iqtree2 -redo --polytomy -bb {int(n_bootstrap)} -m GTR --prefix {prefix} -s {aln_path}'
    else:
        iqtree_cmd = iqtree_cmd = f'iqtree2 -redo --polytomy -m GTR --prefix {prefix} -s {aln_path}'
    subprocess.run(shlex.split(iqtree_cmd))
    return(prefix+'.treefile')
'''


def align_time_tree(aln_path, tre_path, date_dict, root_name, resolve_polytomies, ci_width):
    seqs = list(SeqIO.parse(aln_path, 'fasta'))
    # remove spaces in names for IQtree
    tre = Phylo.read(tre_path, 'newick')
    # dates are for ClockTree class, all other ClockTree parameters are default
    # tree, aln, gtr arguments are passed to TreeAnc class
    tt = TreeTime(tree=tre, aln=aln_path, ref=None, dates=date_dict,
                  gtr='JC69', precision='auto')
    # treetime clock filter will mark, but not remove bad tips
    if root_name != 'best':
        root_names = [i.name for i in tt.tree.get_terminals() if root_name in i.name]
        if len(root_names) > 1:
            raise Exception('multiple root matches')
        root_name = root_names[0]
    tt.clock_filter(reroot=root_name, n_iqd=4, plot=False) 
    # remove them explicitly
    leaves = [x for x in tt.tree.get_terminals()]
    bad_seqs = [item.name for item in leaves if item.bad_branch==True]
    print(f'{len(bad_seqs)} tips failed clock filter')
    pd.DataFrame(bad_seqs).to_csv(tre_path.split('.')[0]+'_bad_seqs.tsv', sep='\t', header=None, index=None)
    seqs = [item for item in seqs if item.id not in bad_seqs]
    time_filtered_seqs_path = tre_path.split('.')[0]+'_time.fasta'
    with open(time_filtered_seqs_path, 'w') as out:
        SeqIO.write(seqs, out, 'fasta')
    # now remove them from tree
    remove_count = 0
    for n in leaves:
        if n.bad_branch:
            remove_count += 1
            tt.tree.prune(n)
            print('pruning leaf ', n.name)
    print(f'{remove_count} sequences removed because they\'re too divergent')
    tt.prepare_tree()
    print(resolve_polytomies)
    tt.run(root=root_name, Tc='skyline', time_marginal='assign', 
        resolve_polytomies=resolve_polytomies, 
        max_iter=4, fixed_clock_rate=0.001, 
        vary_rate=0.0005, use_covariation=True)
    '''
    # modified
    tt.run(root=root_name, Tc='skyline', time_marginal=True, 
        resolve_polytomies=True, 
        max_iter=2, fixed_clock_rate=0.001, 
        vary_rate=0.0005, use_covariation=True)
        
    # original
    tt.run(infer_gtr=True, root=root_name, Tc='skyline', time_marginal='assign', 
        branch_length_mode='auto', resolve_polytomies=True, 
        max_iter=2, fixed_pi=None, fixed_clock_rate=0.001, 
        vary_rate=0.0005, use_covariation=True)
    '''
    times = pd.DataFrame({'name': [item.name for item in tt.tree.find_clades()],
                          'date': [item.numdate for item in tt.tree.find_clades()],
                          f'lower ({ci_width})': 
                            [list(tt.get_max_posterior_region(item, ci_width))[0]
                                for item in tt.tree.find_clades()],
                          f'upper ({ci_width})': 
                            [list(tt.get_max_posterior_region(item, ci_width))[1]
                                for item in tt.tree.find_clades()]}, 
                          index = range(0, len([item for item in tt.tree.find_clades()])))
    times.to_csv(tre_path.split('.')[0]+'_refined_node_times.csv', index=False)
    for n in tt.tree.find_clades():
        n.branch_length = n.mutation_length
    refined_mut_tree_path = tre_path.split('.')[0]+'_refine.nwk'
    with open(refined_mut_tree_path, 'w') as out:
        Phylo.write(tt.tree, out, 'newick', format_branch_length="%1.5f")
    tt.branch_length_to_years()
    refined_tree_path = tre_path.split('.')[0]+'_refine_time.nwk'
    with open(refined_tree_path, 'w') as out:
        Phylo.write(tt.tree, out, 'newick', format_branch_length="%1.5f")
    return(time_filtered_seqs_path, tt)


# this is sloppy, improve
def make_date_dict(seqs_path, seq_name_sep, seq_name_field, metadata_path, 
  metadata_delim, metadata_name_col, metadata_date_col, metadata_date_fmt):
    from utils import numeric_from_datetime, convert_to_datetime
    metadata=pd.read_csv(metadata_path, 
        sep=metadata_delim, header=None)
    metadata[metadata_date_col] = \
        metadata[metadata_date_col].apply(lambda k: convert_to_datetime(k, metadata_date_fmt))
    metadata['numeric_date'] = \
        metadata.loc[metadata[metadata_date_col].notnull(), metadata_date_col].apply(numeric_from_datetime)
    seqs = list(SeqIO.parse(seqs_path, 'fasta'))
    seq_name_dict = \
        {i.description.split(seq_name_sep)[seq_name_field]: i.description for i in seqs}
    metadata['seq_name'] = \
        metadata[metadata_name_col].apply(lambda k: seq_name_dict[k] if k in seq_name_dict.keys() else np.nan)
    date_dict = {i['seq_name']: i['numeric_date'] for idx, i in metadata.iterrows()}
    return(date_dict)


def run():
    parser = argparse.ArgumentParser()
    # input files
    parser.add_argument('--outdir', default='results')
    parser.add_argument('--sequences', 
        default=None)
    parser.add_argument('--seqNameSep',
        default='|')
    parser.add_argument('--seqNameField',
        default=1, type=int, 
        help='which field in the sequence name corresponds to the metadata name column')
    parser.add_argument('--metadata', 
        default=None,
        help='metadata file path, assumes no header')
    parser.add_argument('--metadataSep',
        default=',', 
        help='metadata delimiter')
    parser.add_argument('--metadataNameCol',
        default=0,
        type=int,
        help='which column in the metadata file contains sequence names (corresponds to seqNameField)')
    parser.add_argument('--metadataDateCol',
        default=1,
        type=int,
        help='which column in the metadata file contains sequece dates')
    parser.add_argument('--metadataDateFmt',
        default='%Y-%m-%d',
        help='format of the date column')
    parser.add_argument('--tree', 
        default=None)
    parser.add_argument('--resolve_polytomies',
        default='False', type=str)
    parser.add_argument('--root_name', 
        default='best')
    parser.add_argument('--ci_width', 
        default=0.9,
        help='width of the treetime confidence interval',
        type=float)
    args = parser.parse_args()
    if args.resolve_polytomies.lower() == 'false':
        args.resolve_polytomies = False
    else:
        args.resolve_polytomies = True
    #args.sequences = 'data/france_random_global_ref_aligned.fasta'
    #args.tree = 'data/france_random_global_ref_aligned.fasta.treefile'
    #args.metadata = 'data/france_random_global_ref.tsv'
    #args.root_name = 'EPI_ISL_402125'
    #args.metadataNameCol = 1
    #args.metadataDateCol = 2
    #args.metadataSep = '\t'
    date_dict = make_date_dict(args.sequences, args.seqNameSep, 
        args.seqNameField, args.metadata, args.metadataSep, 
        args.metadataNameCol, args.metadataDateCol, args.metadataDateFmt)
    time_filtered_seqs_path, tt = \
           align_time_tree(args.sequences, args.tree, date_dict, 
            args.root_name, args.resolve_polytomies, args.ci_width)  


if __name__ == "__main__":
    run()

