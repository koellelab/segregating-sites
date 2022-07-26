import argparse
from scipy.special import factorial
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np


def poisson(k, lamb):
    """poisson pdf, parameter lamb is the fit parameter"""
    return (lamb**k/factorial(k)) * np.exp(-lamb)


def negative_log_likelihood(params, data):
    """
    The negative log-Likelihood-Function
    """
    lnl = - np.sum(np.log(poisson(data, params[0])))
    return lnl


def poissfit(dists, a=0.05):
	from scipy.stats import chi2
	mu = dists.mean()
	n_events = dists.sum()
	n_trials = dists.size
	lb = (chi2.ppf(a/2, 2*n_events)/2)/n_trials
	ub = (chi2.ppf(1-a/2, 2*(n_events+1))/2)/n_trials
	return(mu, (lb, ub))


def run():
    parser = argparse.ArgumentParser()
    # input files
    parser.add_argument('--distDat', default='combined_metadata_dist.tsv')
    parser.add_argument('--distCol', default=3, type=int)
    args = parser.parse_args()
    #args.distDat = 'data/combined_metadata_dist.tsv'
    dists = pd.read_csv(args.distDat, header=None, sep='\t').values[:,args.distCol].astype(int)
    mu, ci = poissfit(dists)
    print(f'the MLE of lambda is {mu} with a 95% CI between {ci} (exact methd)')
    with open(args.distDat.replace('.tsv', '_mle.tsv'), 'w') as fp:
    	fp.write(f'{mu} {ci}')


    x_s = np.arange(0, dists.max()+2)
    dat_probs = np.bincount(dists, minlength=x_s.size)
    dat_probs = dat_probs/dat_probs.sum()
    fig, ax = plt.subplots(constrained_layout=True)
    ax.bar(x_s, dat_probs, color="#5E81AC", edgecolor='#333333', label=f'data (n={dists.size})', zorder=1)
    ax.scatter(x_s, poisson(x_s, mu), color='#333333', label='ml estimate', zorder=3)
    for idx, x in enumerate(x_s):
        if idx == 0:
            ax.plot((x,x), (poisson(x, ci[0]), poisson(x, ci[1])), color='#333333', zorder=2, lw=2, label='95% CI')
        else:
            ax.plot((x,x), (poisson(x, ci[0]), poisson(x, ci[1])), color='#333333', zorder=2, lw=2)


    ax.legend()
    ax.set_xlabel(r"# of mutations per transmission")
    ax.set_ylabel('probability')
    ax.set_ylim(-0.039, 1.039)
    ax.set_xlim(-1.039, x_s.max()+1+0.039)
    ax.set_xticks(np.arange(0, x_s.max()+1, 1))

    ax.text(0.05, 0.9, f'$\mu$ = {round(mu,2)} {tuple([round(i, 2) for i in ci])}', transform=ax.transAxes)


    fig.savefig('mu_estimates.pdf')
    plt.close()



if __name__ == "__main__":
    run()


