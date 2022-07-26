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


def plot_data(mu_dat, mu_est, out_name):
	fig = plt.figure( 
		figsize=(6.4*2, 4.8*2), 
		constrained_layout=True)
	fig, ax3 = plt.subplots(1,1, figsize=(6.4, 4.8), constrained_layout=True)
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
	
	fig.savefig(f'{out_name}.pdf')
	plt.close()


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--muDat', default=None)
	parser.add_argument('--muEst', default=None)
	parser.add_argument('--outName', default='figure')
	args = parser.parse_args()
	#args.muEst = 'data/combined_metadata_dist_mle.tsv'
	#args.muDat = 'data/combined_metadata_dist.tsv'
	# get_dates
	mu_est = \
		pd.read_csv(args.muEst, header=None, sep=' ')
	mu_est[0] = mu_est[0].astype(float)
	mu_est[1] = mu_est[1].str.replace(',', '').str.replace('(', '').astype(float)
	mu_est[2] = mu_est[2].str.replace(')', '').astype(float)
	mu_est = mu_est.iloc[0,:].values
	mu_dat = pd.read_csv(args.muDat, sep='\t', header= None)
	plot_data(mu_dat, mu_est, args.outName)


if __name__ == "__main__":
    run()

