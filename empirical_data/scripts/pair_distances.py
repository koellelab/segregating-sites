import numpy as np
import pandas as pd
import argparse
import json



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


def format_seqs_arr(s, n_seqs, nuc_dict={97:[97], 99:[99], 116:[116], 103:[103], 110: [97, 99, 116, 103, 110]}):
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


def import_fasta(fasta_path, nuc_dict={97:[97], 99:[99], 116:[116], 103:[103], 110: [97, 99, 116, 103, 110]}):
	s_names = []
	all_s = ''
	fh = open(fasta_path, 'rt')
	with fh as fasta:
	    for h,s in read_fasta(fasta):
	        s_names.append(h)
	        all_s += s
	fh.close()
	s_arr = format_seqs_arr(all_s, len(s_names), nuc_dict=nuc_dict)
	return(s_arr, s_names)


def make_nuc_arr(nuc_dict={97:[97], 99:[99], 116:[116], 103:[103], 110: [97, 99, 116, 103, 110]}):
	uniq_nucs = np.unique(np.hstack([list(nuc_dict.keys()), np.hstack(list(nuc_dict.values()))]))
	uniq_nucs_map = {i: idx for idx, i in enumerate(uniq_nucs)}
	# probably more efficient way to do this
	nuc_arr =  np.zeros((len(uniq_nucs), len(uniq_nucs)))
	for key, value in nuc_dict.items():
	    nuc_arr[uniq_nucs_map[key], [uniq_nucs_map[i] for i in value]] = 1
	    nuc_arr[[uniq_nucs_map[i] for i in value], [uniq_nucs_map[key]]] = 1
	return(nuc_arr, uniq_nucs_map)


def calc_pair_dists(s_arr, s_names, pairs, nuc_dict):
	s_names_dict = {i:idx for idx, i in enumerate(s_names)}
	nuc_arr, uniq_nucs_map = make_nuc_arr(nuc_dict)
	# todo confirm that symmetric is correct
	u,inv = np.unique(s_arr,return_inverse = True)
	mapped_seqs = np.array([uniq_nucs_map[x] for x in u])[inv].reshape(s_arr.shape)
	# may be possible to do this faster without a loop?
	# todo improve
	dists = np.zeros(len(pairs), dtype='int32')
	for idx,i in enumerate(pairs):
		a = s_names_dict[i[0]]
		b = s_names_dict[i[1]]
		dists[idx] = np.where(nuc_arr[mapped_seqs[a,:],mapped_seqs[b,:]] == 0)[0].shape[0]
	return(dists)


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--seqs')
	parser.add_argument('--pairDat')
	parser.add_argument('--nucDict')
	args = parser.parse_args()
	#args.seqs = 'data/gisaid_hcov-19_2021_06_22_20_ref_aligned_ref_filtered_masked.fasta'
	#args.pairDat = 'data/combined_metadata.tsv'
	#args.nucDict = 'scripts/nuc_dict_all.json'
	pair_dat = pd.read_csv(args.pairDat, sep='\t', header=None).dropna()
	nuc_dict = \
		    {np.frombuffer(key.lower().encode(), dtype=np.int8)[0]: 
		    	np.frombuffer(''.join(value).lower().encode(), dtype=np.int8) 
		    	for key, value in json.load(open(args.nucDict, 'r')).items()}
	s_arr, s_names = import_fasta(args.seqs, nuc_dict=nuc_dict)
	s_names = [i.split('|')[1] for i in s_names]

	pair_dat['dist'] = \
		calc_pair_dists(s_arr, s_names, 
			[tuple(i[1:]) for i in pair_dat.values], nuc_dict)

	pair_dat.to_csv(args.pairDat.replace('.tsv', '_dist.tsv'), header=None, index=None, sep='\t')



if __name__ == "__main__":
    run()



