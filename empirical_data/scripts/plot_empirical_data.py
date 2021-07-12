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


def poisson(k, lamb):
    """poisson pdf, parameter lamb is the fit parameter"""
    return (lamb**k/factorial(k)) * np.exp(-lamb)


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


def get_date_dict(metadata_path):
	metadata = pd.read_csv(metadata_path, header=None, sep='\t')
	metadata[2] = pd.to_datetime(metadata[2])
	date_dict = {i[1]: i[2] for i in metadata.values}
	return(date_dict)


def get_dist_dat(dist_path, date_dict):
	dist_dat = pd.read_csv(dist_path, header=None, sep='\t')
	dist_dat[dist_dat.shape[1]] = dist_dat[0].apply(lambda k: date_dict[k.split('|')[1]])
	# regress on the dist dat
	regress = linregress(dist_dat[dist_dat.shape[1]-1].apply(numeric_from_datetime), 
		dist_dat[1])
	return(dist_dat, regress)


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


def plot_data(dist_dat, regress, s_dat, mu_dat, mu_est, out_name):
	regress_range = \
			(dist_dat[dist_dat.shape[1]-1].min(), 
				dist_dat[dist_dat.shape[1]-1].max())
	fig = plt.figure( 
		figsize=(6.4*2, 4.8*2), 
		constrained_layout=True)
	fig, axs = plt.subplots(2,2, figsize=(6.4*2, 4.8*2), constrained_layout=True)
	ax0 = axs[0,0]
	ax1 = axs[1, 0]
	ax2 = axs[0, 1]
	ax0, ax1, ax2, ax3 = axs.flatten()
	ax0.scatter(dist_dat[2], dist_dat[1], 
		alpha=0.25, color='#81a1c1')
	ax0.plot(regress_range, 
		[regress.slope*(numeric_from_datetime(x)) + 
			regress.intercept for x in regress_range], 
		color='#BF616A')
	ax0.text(0.05, 0.8, 
		str(round(regress.slope)) + " subs/yr, \n{:.2e} subs/site/yr".format(regress.slope/29903),
		transform=ax0.transAxes, color='#BF616A', size=12)
	ax0.set_ylabel('distance to Wuhan/Hu-1')
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
	# mu dat
	x_s = np.arange(0, mu_dat[3].max()+2)
	dist_counts = np.bincount(mu_dat[3].values,  minlength=x_s.size)
	dist_probs = dist_counts/dist_counts.sum()
	ax3.bar(x_s, dist_probs, color="#5E81AC", edgecolor='#333333', label=f'data (n={mu_dat[3].size})', zorder=1)
	mu = mu_est[0]
	ci = mu_est[1:]
	ax3.scatter(x_s, poisson(x_s, mu), color='#333333', label='ml estimate', zorder=3)
	for idx, x in enumerate(x_s):
		if idx == 0:
			ax3.plot((x,x), (poisson(x, ci[0]), poisson(x, ci[1])), color='#333333', zorder=2, lw=2, label='95% CI')
		else:
			ax3.plot((x,x), (poisson(x, ci[0]), poisson(x, ci[1])), color='#333333', zorder=2, lw=2)
	ax3.legend()
	ax3.set_xlabel(r"# of mutations per transmission")
	ax3.set_ylabel('probability')
	ax3.set_ylim(-0.039, 1.039)
	ax3.set_xlim(-1.039, x_s.max()+1+0.039)
	ax3.set_xticks(np.arange(0, x_s.max()+1, 1))
	ax3.text(0.05, 0.9, f'$\mu$ = {round(mu,2)} {tuple([round(i, 2) for i in ci])}', transform=ax3.transAxes)
	time_axis = [ax0, ax1, ax2]
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
	for ax in time_axis[:-1]:
		ax.set_xlabel(None)
		ax.set_xticklabels([])
	for ax_idx, ax in enumerate([ax0, ax1, ax2]):
		ax.text(-0.1, 1.05, string.ascii_uppercase[ax_idx], transform=ax.transAxes, 
			size=16, weight='bold', va="top")
	fig.savefig(f'{out_name}.pdf')
	plt.close()


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--distDat', default=None)
	parser.add_argument('--sDat', default=None)
	parser.add_argument('--muDat', default=None)
	parser.add_argument('--muEst', default=None)
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
	date_dict = get_date_dict(args.metadata)
	dist_dat, regress = \
		get_dist_dat(args.distDat, date_dict)
	s_dat = get_s_dat(args.sDat)
	mu_est = \
		pd.read_csv(args.muEst, header=None, sep=' ')
	mu_est[0] = mu_est[0].astype(float)
	mu_est[1] = mu_est[1].str.replace(',', '').str.replace('(', '').astype(float)
	mu_est[2] = mu_est[2].str.replace(')', '').astype(float)
	mu_est = mu_est.iloc[0,:].values
	mu_dat = pd.read_csv(args.muDat, sep='\t', header= None)
	plot_data(dist_dat, regress, s_dat, mu_dat, mu_est, args.outName)


if __name__ == "__main__":
    run()

