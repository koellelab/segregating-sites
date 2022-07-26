import argparse
import pandas as pd
import numpy as np
import ast
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



def make_nuc_arr(nuc_dict={97:[97], 99:[99], 116:[116], 103:[103], 110: [97, 99, 116, 103, 110]}):
	uniq_nucs = np.unique(np.hstack([list(nuc_dict.keys()), np.hstack(list(nuc_dict.values()))]))
	uniq_nucs_map = {i: idx for idx, i in enumerate(uniq_nucs)}
	# probably more efficient way to do this
	nuc_arr =  np.zeros((len(uniq_nucs), len(uniq_nucs)))
	for key, value in nuc_dict.items():
	    nuc_arr[uniq_nucs_map[key], [uniq_nucs_map[i] for i in value]] = 1
	    #nuc_arr[[uniq_nucs_map[i] for i in value], [uniq_nucs_map[key]]] = 1
	return(nuc_arr, uniq_nucs_map)




def map_arr(arr, map_dict):
    k = np.array(list(map_dict.keys()))
    v = np.array(list(map_dict.values()))
    mapping_arr = np.zeros(k.max()+1,dtype=v.dtype) #k,v from approach #1
    mapping_arr[k] = v
    mapped_arr = mapping_arr[arr]
    return(mapped_arr)



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


def calc_seg_sites(seqs, nuc_dict={97:[97], 99:[99], 116:[116], 103:[103], 110: [97, 99, 116, 103, 110]}):
	if seqs.shape[0] >= 1:
		u,inv = np.unique(seqs,return_inverse = True)
		nucs = np.array([nuc_dict[x] for x in u], dtype=object)[inv].reshape(seqs.shape)
		l = max(max(nuc_dict.keys()), max([max(item) for item in nuc_dict.values()]))
		bins = np.array([np.bincount(np.concatenate(i), minlength=l).max() for i in nucs.T])
		seg_sites = (np.where(bins < seqs.shape[0])[0])
		return(seg_sites)
	else:
		return(np.nan)


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--toMatch', type=int)
	parser.add_argument('--maxDist', type=int)
	parser.add_argument('--ccSNPs')
	parser.add_argument('--cc')
	parser.add_argument('--seqs')
	parser.add_argument('--dists')
	parser.add_argument('--nucDict')
	parser.add_argument('--ref')
	args = parser.parse_args()
	#args.ccSNPs = 'data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_snp.tsv'
	#args.cc = 'data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc.tsv'
	#args.dists = 'data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref.dist'
	#args.nucDict = 'scripts/nuc_dict_all.json'
	#args.ref = 'data/EPI_ISL_402125.fasta'
	#args.toMatch= 0
	#args.maxDist = 7
	# read in nucleotide dictionary to handle ambiguities
	nuc_dict = \
	    {np.frombuffer(key.lower().encode(), dtype=np.int8)[0]: 
	    	np.frombuffer(''.join(value).lower().encode(), dtype=np.int8) 
	    	for key, value in json.load(open(args.nucDict, 'r')).items()}
	# read in reference sequence
	ref_seq_arr = import_fasta(args.ref)[0][0]
	# read in SNPs that define connected components
	cc_snps = pd.read_csv(args.ccSNPs, header=None, sep='\t')
	cc_snps_dict = {i[0]: [ast.literal_eval(k) for k in set(i[1:])-set([np.nan])] for i in cc_snps.values}
	# generate snp alignment from connected component snps
	snp_aln = np.tile(ref_seq_arr, cc_snps.shape[0]).reshape((cc_snps.shape[0], ref_seq_arr.size))
	for key_idx, key in enumerate(cc_snps_dict):
		for snp in cc_snps_dict[key]:
			snp_aln[key_idx, int(snp[1])-1] = np.frombuffer(snp[2].lower().encode(), dtype=np.int8)
	# get sites that are polymorphic in the requested cc id 
	to_match_var_sites = calc_seg_sites(np.vstack([snp_aln[args.toMatch,:], ref_seq_arr]), nuc_dict=nuc_dict)
	# alignment of just those sites
	to_match_snp_aln = snp_aln[:,to_match_var_sites]
	# now get any other CCs which match at these sites
	match = np.where(calc_pairwise_dists(to_match_snp_aln, nuc_dict=nuc_dict)[args.toMatch,:] == 0 )[0]
	# read in CC assignments
	cc = pd.read_csv(args.cc, header=None, sep='\t')
	match_seqs = cc[cc[1].isin(match)][0].values
	# read in pairwise distances
	dists = pd.read_csv(args.dists, sep='\t', header=None)
	dists = dists.set_index(0, drop=True)
	keep_ids = \
		[idx for idx, i in enumerate(dists.index) if i in match_seqs]
	dists = dists.iloc[keep_ids, keep_ids]
	dists.values[np.diag_indices(dists.shape[0])] = 10000000
	match_seqs = [i for i in match_seqs if dists.loc[i,:].min() <= args.maxDist]
	# save to file
	with open(args.cc.replace('.tsv', '_match_largest.tsv'), 'w') as fp:
		for i in match_seqs:
			fp.write(i+'\n')



if __name__ == "__main__":
    run()



