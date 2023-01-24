import numpy as np
import argparse
import json
import pandas as pd
import ast
import matplotlib.pyplot as plt
from collections import Counter

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))



def import_fasta(fasta_path, nuc_dict={97:[97, 110], 99:[99, 110], 116:[116, 110], 103:[103, 110]}):
	s_names = []
	all_s = ''
	fh = open(fasta_path, 'rt')
	with fh as fasta:
	    for h,s in read_fasta(fasta):
	        s_names.append(h)
	        all_s += s
	fh.close()
	s_arr = format_seqs_arr(all_s, len(s_names), nuc_dict=nuc_dict)
	return(s_arr, np.array(s_names))



def format_seqs_arr(s, n_seqs, nuc_dict={97:[97, 110], 99:[99, 110], 116:[116, 110], 103:[103, 110]}):
    n_seqs = int(n_seqs)
    size = int(len(s)/n_seqs)
    seqs_arr = \
        np.frombuffer(s.lower().encode(), dtype=np.int8)
    seqs_arr = np.copy(seqs_arr)
    #[(97, 'A'), (114, 'R'), (119, 'W'), (109, 'M'), (100, 'D'), (104, 'H'), (118, 'V'), 
    #(110, 'N'), (99, 'C'), (121, 'Y'), (115, 'S'), (109, 'M'), (98, 'B'), (104, 'H'), 
    #(118, 'V'), (110, 'N'), (117, 'U'), (121, 'Y'), (119, 'W'), (107, 'K'), (98, 'B'), 
    #(100, 'D'), (104, 'H'), (110, 'N'), (103, 'G'), (114, 'R'), (115, 'S'), (107, 'K'), 
    #(98, 'B'), (100, 'D'), (118, 'V'), (110, 'N'), (116, 'T')]
    all_nucs = np.hstack(list(nuc_dict.values()))
    seqs_arr[~np.in1d(seqs_arr, all_nucs)] = 110
    seqs_arr = \
        seqs_arr.reshape((n_seqs, int(seqs_arr.shape[0]/n_seqs)))
    return(seqs_arr)


def calc_pairwise_dists(seqs, nuc_dict={97:[97, 110], 99:[99, 110], 116:[116, 110], 103:[103, 110]}):	
	uniq_nucs = sorted(list(set([item for sublist in nuc_dict.values() for item in sublist])))
	uniq_nucs_map = {i: idx for idx, i in enumerate(uniq_nucs)}
	nuc_arr =  np.zeros((len(uniq_nucs), len(uniq_nucs)))
	for key, value in nuc_dict.items():
		nuc_arr[uniq_nucs_map[key], [uniq_nucs_map[i] for i in value]] = 1
		nuc_arr[[uniq_nucs_map[i] for i in value], [uniq_nucs_map[key]]] = 1
	u,inv = np.unique(seqs,return_inverse = True)
	mapped_seqs = np.array([uniq_nucs_map[x] for x in u])[inv].reshape(seqs.shape)
	# may be possible to do this faster without a loop?
	# todo improve
	combos = list(itertools.combinations(np.arange(seqs.shape[0]), 2))
	dist_arr = np.zeros((seqs.shape[0], seqs.shape[0]))
	for a,b in combos:
		d = np.where(nuc_arr[mapped_seqs[a,:],mapped_seqs[b,:]] == 0)[0].shape[0]
		dist_arr[a,b] = d
		dist_arr[b,a] = d
	return(dist_arr)


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--SNPs')
	parser.add_argument('--seqs')
	parser.add_argument('--nucDict')
	parser.add_argument('--ref')
	parser.add_argument('--outName')
	args = parser.parse_args()
	#args.SNPs = 'test.tsv'
	#args.nucDict = 'scripts/nuc_dict_all.json'
	#args.ref = 'data/EPI_ISL_402125.fasta'
	#args.seqs = 'data/gisaid_hcov-19_2021_06_12_01_aligned_ref_filtered_masked.fasta'
	#args.outName = '_match_largest_snps'
	nuc_dict = \
	    {np.frombuffer(key.lower().encode(), dtype=np.int8)[0]: 
	    	np.frombuffer(''.join(value).lower().encode(), dtype=np.int8) 
	    	for key, value in json.load(open(args.nucDict, 'r')).items()}
	# read in reference sequence
	ref_seq_arr = import_fasta(args.ref, nuc_dict=nuc_dict)[0][0]
	# read in SNPs that define connected components
	snps = pd.read_csv(args.SNPs, header=None, sep='\t').T.dropna()
	snps = snps[0].apply(ast.literal_eval)
	snps_locs = \
		np.array([i[1]-1 for i in snps if i[2].lower() != 'n'])
	snps_nucs = \
		np.frombuffer(''.join([i[2]for i in snps if i[2].lower() != 'n']).lower().encode(), dtype=np.int8)
	mask = np.ones(ref_seq_arr.shape[0], dtype='bool')
	mask[snps_locs] = False
	# now read in exogeneous sequences
	seqs_arr, seqs_names = import_fasta(args.seqs, nuc_dict=nuc_dict)
	# get sequences which match at the clade defining sites
	#match_seqs_names = \
	#	pd.DataFrame(seqs_names[((seqs_arr[:,snps_locs] == snps_nucs).all(axis=1) & \
	#	((seqs_arr[:,mask] == ref_seq_arr[mask]) | \
	#		(seqs_arr[:,mask] == 110)).all(axis=1))])
	match_seqs_names = \
		pd.DataFrame(seqs_names[(seqs_arr[:,snps_locs] == snps_nucs).all(axis=1)])

	match_seqs_names['date'] = \
		pd.to_datetime(match_seqs_names[0].str.split('|', expand=True).iloc[:,-1])

	match_seqs_names.sort_values(by='date').to_csv('.'.join(args.seqs.split('.')[:-1])+'_'+args.outName+'.tsv', sep='\t', header=None)

	counted_dates = Counter(match_seqs_names['date'])

	import matplotlib.dates as mdates
	myFmt = mdates.DateFormatter('%b. %d')



	fig, ax = plt.subplots()
	ax.bar(counted_dates.keys(), counted_dates.values(), facecolor='steelblue', edgecolor='#333333')
	ax.set_xlabel('date (2020)')
	ax.set_ylabel('count')
	ax.xaxis.set_major_formatter(myFmt)
	fig.savefig('basal_sampling_dates.pdf')



if __name__ == "__main__":
	run()

