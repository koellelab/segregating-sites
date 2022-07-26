import scipy.sparse
import pandas as pd
import numpy as np
import argparse
import copy
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



def import_sparse_arr(file_path, snp_threshold, skip_ids=None):
	if type(skip_ids) != type(None):
		skip_ids_arr = np.array(skip_ids)
		skip_ids = set(skip_ids)
	else:
		skip_ids_arr = np.array([])
		skip_ids = set([])
	names = []
	rows = []
	cols = []
	line_len = 0
	with open(file_path) as f:
		for line_idx, line in enumerate(f):
			if line_idx not in skip_ids:
				split_line = line.strip().split()
				names.append(split_line[0])
				line_arr = np.array(split_line[1:], dtype=np.int16)
				connected_indices = np.where(line_arr <= snp_threshold)[0]
				connected_indices = connected_indices[~np.in1d(connected_indices, skip_ids_arr)]
				line_len = line_arr.shape[0]
				cols.extend(connected_indices)
				rows.extend([line_idx]*len(connected_indices))
	dat = np.ones(shape=len(cols), dtype=np.int8)
	sparse_mat = scipy.sparse.csr_matrix((dat, (rows, cols)), 
		shape=(line_len, line_len))
	return(sparse_mat, names)	


def arr_to_fasta(arr):
	nuc_dict = {110: 'N', 97: 'A', 99: 'C', 103: 'G', 116: 'T'}
	keys, inv = np.unique(arr, return_inverse=True)
	vals = np.array([nuc_dict[key] for key in keys])
	nucs = vals[inv]
	return(nucs.reshape(arr.shape))


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


# this is slow on big alignments, improve
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


def get_shared_snps(aln, ref=None, nuc_dict={97:[97, 110], 99:[99, 110], 116:[116, 110], 103:[103, 110]}):	
	# returns 1 indexed SNP positions
	# get segregating sites amongst the aln
	aln_segs = calc_seg_sites(aln, nuc_dict=nuc_dict)
	# get segregating sites amongst the aln + ref
	aln_ref_segs = calc_seg_sites(np.vstack([aln, ref]), nuc_dict=nuc_dict)
	# the difference between these two are the shared SNPs
	shared_snps = np.setdiff1d(aln_ref_segs, aln_segs)
	# most common nucleotide at each shared SNP site
	# Counter returns sorted
	from collections import Counter
	shared_snps_alt = [list(Counter(aln[:,i]).keys())[0] for i in shared_snps]
	# shared_snps any sites where everything in aln is N
	all_N = np.where((aln == 110).all(axis=0))[0]
	shared_snps = np.concatenate([shared_snps, all_N])
	shared_snps_alt.extend(['110'] * all_N.shape[0])
	shared_snps_alt = np.array([shared_snps_alt])
	shared_snps_ref = np.array([ref[i] for i in shared_snps])
	sort_idx = np.argsort(shared_snps)
	format_out = list(zip(list(shared_snps_ref.tobytes().decode('UTF-8').upper()),
		shared_snps + 1, 
		list(shared_snps_alt.astype(np.int8).tobytes().decode('UTF-8').upper())))
	return(format_out)


def calc_minimum_distance(arr_1, arr_2):
	snp_dict = {}
	arr_1_dist = np.inf
	for arr_1_idx, arr_1_seq in enumerate(arr_1):
		for arr_2_idx, arr_2_seq in enumerate(arr_2):
			snps = ((arr_1_seq != arr_2_seq) & 
						(arr_1_seq != 110) & 
						(arr_2_seq != 110))
			dist = snps.sum()
			arr_1_dist = min(arr_1_dist, dist)
			snp_dict[(arr_1_idx, arr_2_idx)] = np.where(snps==True)[0]
	return(arr_1_dist, snp_dict)



def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--dist')
	parser.add_argument('--seqs')
	parser.add_argument('--metadata')
	parser.add_argument('--ref')
	parser.add_argument('--snpThreshold', default=1, type=int)
	parser.add_argument('--nucDict', default=None)
	args = parser.parse_args()
	#args.dist = 'data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref.dist'
	#args.metadata = 'data/gisaid_hcov-19_2021_04_29_16.tsv'
	#args.ref = 'data/EPI_ISL_402125.fasta'
	#args.seqs = 'data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref.fasta'
	#args.nucDict = 'scripts/nuc_dict_all.json'
	nuc_dict = \
	    {np.frombuffer(key.lower().encode(), dtype=np.int8)[0]: 
	    	np.frombuffer(''.join(value).lower().encode(), dtype=np.int8) 
	    	for key, value in json.load(open(args.nucDict, 'r')).items()}
	metadata = pd.read_csv(args.metadata, sep='\t')
	metadata['Collection date'] = pd.to_datetime(metadata['Collection date'])
	date_dict = {i[0].replace(' ', ''): i[2] for i in metadata.values}
	# read in reference sequence
	ref_seq_arr = import_fasta(args.ref)[0][0]
	sparse_dist_arr, dist_df_names = import_sparse_arr(args.dist, args.snpThreshold)
	# 1 if connected (<= snp_threshold)
	# 0 if unconnected (> snp_threshold)
	cc = scipy.sparse.csgraph.connected_components(sparse_dist_arr)
	cc_df = pd.DataFrame([[dist_df_names[idx], i] for idx, i in enumerate(cc[1])])
	cc_df.sort_values(by=1).to_csv(args.dist.replace('.dist', '_cc.tsv'), sep='\t', header=None, index=None)
	# read in sequences
	s_arr, s_names = import_fasta(args.seqs, nuc_dict=nuc_dict)
	# converts names to appropriate format
	s_names_format = [i.split('|')[0] for i in s_names]
	# find defining SNPs/N sites for each connected component
	cc_shared = {}
	for cc_idx, cc_dat in cc_df.groupby(1):
		print(cc_idx)
		cc_names = cc_dat[0]
		cc_names_idxs = [s_names.index(i) for i in cc_names]
		cc_seqs = s_arr[cc_names_idxs, :]
		if cc_seqs.shape[0] != cc_names.shape[0]:
			raise Exception('not all names found in alignment')
		cc_shared[cc_idx] = get_shared_snps(cc_seqs, ref=ref_seq_arr, nuc_dict=nuc_dict)
	pd.DataFrame.from_dict(cc_shared, orient='index').to_csv(args.dist.replace('.dist', '_cc_snp.tsv'), sep='\t', header=None)



if __name__ == "__main__":
    run()

