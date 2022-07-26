import matplotlib.pyplot as plt
import pandas as pd
import argparse
from scipy.stats import linregress
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
import datetime
import string
import numpy as np
from scipy.special import factorial


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



def get_s_dat(s_path):
	# read in and format s dat
	s_dat = pd.read_csv(s_path, sep='\t')
	s_dat['dt'] = \
		s_dat['window_end'].apply(lambda k: 
			datetime.datetime.fromordinal(int(k-366)))
	# remove blank first line, if present
	if s_dat['n'].iloc[0] == 0:
		s_dat = s_dat.iloc[1:,:]
	# remove blank trailing line, if present
	if s_dat['n'].iloc[-1] == 0:
		s_dat = s_dat.iloc[:-1,:]
	return(s_dat)


def plot_data(s_dat, out_name):
	fig = plt.figure( 
		figsize=(6.4*2, 4.8*2), 
		constrained_layout=True)
	fig, axs = plt.subplots(1,2, figsize=(6.4*2, 4.8), constrained_layout=True)
	ax1, ax2 = axs
	# N data
	ax1.scatter(s_dat['dt'], s_dat['n'], 
		color='#81a1c1')
	ax1.plot(s_dat['dt'], s_dat['n'], 
		color='#81a1c1')
	ax1.set_ylabel('# of sequences')
	# S data
	ax2.scatter(s_dat['dt'], s_dat['s'], 
		color='#81a1c1')
	ax2.plot(s_dat['dt'], s_dat['s'], 
		color='#81a1c1')
	ax2.set_ylabel('segregating sites')
	time_axis = [ax1, ax2]
	min_x_val = min([i.get_xlim()[0] for i in time_axis])
	max_x_val = max([i.get_xlim()[1] for i in time_axis])
	for ax in time_axis:
		ax.set_xlabel('date')
		ax.set_xlim(min_x_val, max_x_val)
		ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
		ticks = ax.get_xticks()
		# make min and max spacing the same
		spacing = max(ax.get_xticks()[0] - min_x_val, 
			max_x_val - ax.get_xticks()[-1])
		ax.set_xlim(ax.get_xticks()[0] - spacing, ax.get_xticks()[-1]+spacing)
		ax.set_xticks(s_dat['dt'])
		_ = [tick.set_rotation(45) for tick in ax.get_xticklabels()]
		_ = [tick.set_horizontalalignment('right') for tick in ax.get_xticklabels()]
	#for ax in time_axis[:-1]:
	#	ax.set_xlabel(None)
	#	ax.set_xticklabels([])
	for ax_idx, ax in enumerate([ax1, ax2]):
		ax.text(-0.1, 1.05, string.ascii_uppercase[ax_idx], transform=ax.transAxes, 
			size=16, weight='bold', va="top")
	fig.savefig(f'{out_name}.pdf')
	plt.close()


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--sDat', default=None)
	parser.add_argument('--metadata', default=None)
	parser.add_argument('--outName', default='figure')
	args = parser.parse_args()
	#args.distDat = 'data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_match_largest_EPI_ISL_402125.tsv'
	#args.sDat = 'data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_match_largest_s.tsv'
	#args.metadata = 'data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_match_largest.tsv'
	#args.outName = 'figures/france_s_dat'
	#args.muEst = 'data/combined_metadata_dist_mle.tsv'
	#args.muDat = 'data/combined_metadata_dist.tsv'
	# get_dates
	s_dat = get_s_dat(args.sDat)
	plot_data(s_dat, args.outName)


if __name__ == "__main__":
    run()

