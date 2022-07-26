import pandas as pd
import numpy as np
import argparse
import json
import itertools



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


# maybe useful code in the future
'''


def calc_dists(seqs, target_seq_index, 
  nuc_dict={97:[97], 99:[99], 116:[116], 103:[103], 110: [97, 99, 116, 103, 110]}):
    uniq_nucs = np.unique(np.hstack(list((nuc_dict.values()))))
    uniq_nucs_map = {i: idx for idx, i in enumerate(uniq_nucs)}
    # probably more efficient way to do this
    nuc_arr =  np.zeros((len(uniq_nucs), len(uniq_nucs)))
    for key, value in nuc_dict.items():
        nuc_arr[uniq_nucs_map[key], [uniq_nucs_map[i] for i in value]] = 1
        nuc_arr[[uniq_nucs_map[i] for i in value], [uniq_nucs_map[key]]] = 1
    k = np.array(list(uniq_nucs_map.keys()))
    v = np.array(list(uniq_nucs_map.values()))
    mapping_arr = np.zeros(k.max()+1,dtype=v.dtype) #k,v from approach #1
    mapping_arr[k] = v
    mapped_target = mapped_seqs[target_seq_index]
    dists = (nuc_arr[mapped_seqs, mapped_target] == 0).sum(axis=1)
    return(dists)
'''

def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--seqs', default=None)
	parser.add_argument('--nucDict', default=None)
	args = parser.parse_args()
	#args.seqs = 'data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref.fasta'
	#args.nucDict = 'scripts/nuc_dict_all.json'
	nuc_dict = \
	    {np.frombuffer(key.lower().encode(), dtype=np.int8)[0]: 
	    	np.frombuffer(''.join(value).lower().encode(), dtype=np.int8) 
	    	for key, value in json.load(open(args.nucDict, 'r')).items()}
	s_arr, s_names = import_fasta(args.seqs, nuc_dict=nuc_dict)
	dist_arr = calc_pairwise_dists(s_arr, nuc_dict=nuc_dict)
	with open(args.seqs.replace('.fasta', '.dist'), 'w') as fp:
		for name, row in zip(s_names, dist_arr):
			fp.write(name+'\t'+'\t'.join([str(int(i)) for i in row])+'\n')




if __name__ == "__main__":
    run()