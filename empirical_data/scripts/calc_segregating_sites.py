import pandas as pd
import numpy as np
import argparse
from utils import numeric_from_datetime, datetime_from_numeric, matlab_ordinal_from_numeric
#from scripts.utils import numeric_from_datetime, datetime_from_numeric, matlab_ordinal_from_numeric
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


def calc_seg_sites(seqs, nuc_dict={97:[97], 99:[99], 116:[116], 103:[103], 110: [97, 99, 116, 103, 110]}):
	if seqs.shape[0] >= 1:
		u,inv = np.unique(seqs,return_inverse = True)
		nucs = np.array([nuc_dict[x] for x in u], dtype=object)[inv].reshape(seqs.shape)
		l = max(max(nuc_dict.keys()), max([max(item) for item in nuc_dict.values()]))
		bins = np.array([np.bincount(np.concatenate(i), minlength=l).max() for i in nucs.T])
		seg_sites = (np.where(bins < seqs.shape[0])[0])
		return(seg_sites.shape[0])
	else:
		return(np.nan)



def calc_s_traj(seqs, dates, times, n_per_window=None, 
  nuc_dict={97:[97], 99:[99], 116:[116], 103:[103], 110: [97, 99, 116, 103, 110]}):
	# ASSUMES THAT DATES ARE IN THE SAME ORDER AS THE SEQUENCE ARRAY
	# binning windows this way, each window includes times up to but not including
	# the time at that index. For example if 
	# times = [737812, 737816]
	# then 737811 = window 0, 737812 = window 1, 737814 = window 1, 737816 = window 2
	binned_dates = np.digitize(dates, times)
	# get the indices to sort bins and seqs
	sorted_idxs = np.argsort(binned_dates)
	sorted_bins = binned_dates[sorted_idxs]
	sorted_seqs = seqs[sorted_idxs]
	bin_counts = \
				np.bincount(sorted_bins, 
					minlength=len(times+1))
	split_seqs = \
			np.split(sorted_seqs, 
				np.cumsum(bin_counts))[:-1]
	if n_per_window:
		split_seqs = \
			[i[np.random.default_rng().choice(range(i.shape[0]), 
				size=n_per_window, replace=False)] if i.shape[0] >= n_per_window else np.array([])
				for idx, i in enumerate(split_seqs)]
	bin_s = [calc_seg_sites(i, nuc_dict=nuc_dict) for i in split_seqs]
	return(bin_counts, bin_s)


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--metadata', default=None)
	parser.add_argument('--metadataDelim', default='\t')
	parser.add_argument('--metadataNameCol', default=0, type=int)
	parser.add_argument('--metadataDateCol', default=1, type=int)
	parser.add_argument('--metadataDateFmt', default='%Y-%m-%d')
	parser.add_argument('--tStep', nargs='+', 
	    type=int, default=4)
	parser.add_argument('--windowStart', 
	    type=int, default=-2, help='sequence start date relative to date of first sequence')
	parser.add_argument('--seqs', 
	    default=None)
	parser.add_argument('--seqNameSep', 
	    default='|')
	parser.add_argument('--seqNameField', 
	    default=1, type=int)
	parser.add_argument('--nucDict', 
	    default='scripts/nuc_dict_all.json',
	    help='which nucleotide dictionary to use')
	parser.add_argument('--matlabDate', 
		dest='matlab_date', action='store_true')
	parser.add_argument('--numericDate', 
		dest='matlab_date', action='store_false')
	parser.add_argument('--nPerWindow',  
	    type=int)
	parser.set_defaults(matlab_date=True)
	args = parser.parse_args()
	#args.metadata = 'data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_match_largest.tsv'
	#args.seqs = 'data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_match_largest.fasta'
	#args.metadataDelim = '\t'
	#args.metadataNameCol = 1
	#args.metadataDateCol = 2
	# imorts nucleotide dictionary
	nuc_dict = \
	    {np.frombuffer(key.lower().encode(), dtype=np.int8)[0]: 
	    	np.frombuffer(''.join(value).lower().encode(), dtype=np.int8) 
	    	for key, value in json.load(open(args.nucDict, 'r')).items()}
	# imports fasta alignment
	s_arr, s_names = import_fasta(args.seqs, nuc_dict=nuc_dict)
	# imports date metadata
	print(args.metadataDelim)
	print(args.metadata)
	print(pd.read_csv(args.metadata, sep=args.metadataDelim, header=None))
	dates_df = \
	    pd.read_csv(args.metadata, sep=args.metadataDelim, header=None)\
	        [[args.metadataNameCol, args.metadataDateCol]]
	print(dates_df)
	dates_df[args.metadataDateCol] = \
	    pd.to_datetime(dates_df[args.metadataDateCol], 
	    	format=args.metadataDateFmt).apply(numeric_from_datetime)
	date_dict = {i[0]: i[1] for i in dates_df.values}
	if args.matlab_date:
	    date_dict = {key: int(matlab_ordinal_from_numeric(value)) for key, value in date_dict.items()}
	# creates window times
	dates = [date_dict[i.split(args.seqNameSep)[args.seqNameField]] for i in s_names]
	# times array indicates the first date in each windows
	# thus, we want the last value in the times array to be
	# the maximum observed date + 1 
	# also, numpy doesn't include the last value
	# in arange, so we add a tstep
	# sloppy, fix this 
	print(min(dates))
	print(max(dates))
	n_windows = np.ceil((max(dates)+1-min(dates))/args.tStep)
	t0 = max(dates) + 1 - args.tStep*n_windows
	max_times_value = t0 + n_windows*args.tStep + 1
	times = np.arange(t0, max_times_value, args.tStep)
	#times = np.arange(min(dates) + args.windowStart, max(dates)+args.tStep*2, args.tStep)
	n_traj, s_traj = calc_s_traj(s_arr, dates, times, 
		nuc_dict=nuc_dict,n_per_window=args.nPerWindow)
	# need to adjust times by 1 by how binning occurs
	n_s_traj = pd.DataFrame([times-1, n_traj, [i for i in s_traj]], 
		index=['window_end', 'n', 's']).transpose()
	if max(n_s_traj['window_end']) != max(dates):
		raise Exception('times array mismatch')
	if args.nPerWindow:
		n_s_traj.to_csv(args.seqs.replace('.fasta', f'_s_{args.nPerWindow}.tsv'), sep='\t', index=None)
	else:
		n_s_traj.to_csv(args.seqs.replace('.fasta', f'_s.tsv'), sep='\t', index=None)
	print(list(n_s_traj.columns))
	for idx, row in n_s_traj.iterrows():
		print(row.values)



if __name__ == "__main__":
    run()

