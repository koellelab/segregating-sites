import numpy as np
import argparse
import json
import itertools
import pandas as pd

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
	return(s_arr, s_names)



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


def calc_dists(seqs, ref, nuc_dict={97:[97, 110], 99:[99, 110], 116:[116, 110], 103:[103, 110]}):	
	uniq_nucs = sorted(list(set([item for sublist in nuc_dict.values() for item in sublist])))
	uniq_nucs_map = {i: idx for idx, i in enumerate(uniq_nucs)}
	nuc_arr =  np.zeros((len(uniq_nucs), len(uniq_nucs)))
	for key, value in nuc_dict.items():
		nuc_arr[uniq_nucs_map[key], [uniq_nucs_map[i] for i in value]] = 1
		nuc_arr[[uniq_nucs_map[i] for i in value], [uniq_nucs_map[key]]] = 1
	u,inv = np.unique(seqs,return_inverse = True)
	mapped_seqs = np.array([uniq_nucs_map[x] for x in u])[inv].reshape(seqs.shape)
	u, inv = np.unique(ref, return_inverse= True)
	mapped_ref = np.array([uniq_nucs_map[x] for x in u])[inv].reshape(ref.shape)
	mapped_ref_rep = np.tile(mapped_ref, seqs.shape[0])
	# may be possible to do this faster without a loop?
	# todo improve
	mapped_seqs_flat = mapped_seqs.flatten()
	similarity = nuc_arr[mapped_seqs_flat, mapped_ref_rep].reshape(seqs.shape)
	diffs = (similarity != 1).sum(axis=1)
	return(diffs)


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--querySeqs', default=None)
	parser.add_argument('--targetSeq', default=None)
	parser.add_argument('--nucDict', default=None)
	args = parser.parse_args()
	#args.querySeqs = 'data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_match_largest.fasta'
	#args.targetSeq = 'data/EPI_ISL_402125.fasta'
	#args.nucDict = 'scripts/nuc_dict_all.json'
	nuc_dict = \
	    {np.frombuffer(key.lower().encode(), dtype=np.int8)[0]: 
	    	np.frombuffer(''.join(value).lower().encode(), dtype=np.int8) 
	    	for key, value in json.load(open(args.nucDict, 'r')).items()}
	#args.querySeqs = 'data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref.fasta'
	s_arr, s_names = import_fasta(args.querySeqs, nuc_dict=nuc_dict)
	ref_arr, ref_name = import_fasta(args.targetSeq, nuc_dict=nuc_dict)
	ref_arr = ref_arr[0]
	dist_arr = calc_dists(s_arr, ref_arr, nuc_dict=nuc_dict)
	dist_df = pd.DataFrame([s_names, dist_arr]).T

	target_seq_name = \
		args.targetSeq.split('/')[-1].replace('.fasta', '')
	out_name = \
		args.querySeqs.replace('.fasta', 
			'_' + target_seq_name+'.tsv')
	dist_df.to_csv(out_name, sep='\t', index=None, header=None)



if __name__ == "__main__":
    run()
