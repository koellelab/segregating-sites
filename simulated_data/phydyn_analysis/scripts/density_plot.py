import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import datetime


def datetime_from_numeric(numdate):
    # borrowed from the treetime utilities
    # https://github.com/neherlab/treetime/blob/de6947685fbddc758e36fc4008ddd5f9d696c6d3/treetime/utils.py
    """convert a numeric decimal date to a python datetime object
    Note that this only works for AD dates since the range of datetime objects
    is restricted to year>1.
    Parameters
    ----------
    numdate : float
        numeric date as in 2018.23
    Returns
    -------
    datetime.datetime
        datetime object
    """
    from calendar import isleap
    import datetime
    days_in_year = 366 if isleap(int(numdate)) else 365

    # add a small number of the time elapsed in a year to avoid
    # unexpected behavior for values 1/365, 2/365, etc
    days_elapsed = int(((numdate%1)+1e-10)*days_in_year)
    date = datetime.datetime(int(numdate),1,1) + datetime.timedelta(days=days_elapsed)
    return date


# mostly from arviz
# but works directly with a list
# and also returns median
# https://github.com/arviz-devs/arviz/blob/master/arviz/stats/stats.py
def hpd(dat, qs=[0.025, 0.5, 0.975]):
    width = qs[2] - qs[0]
    # sorts from smallest to largest
    dat = sorted(dat)
    # length of data
    n = len(dat)
    # number of values we are keeping
    # thus, the index of the begining of 
    # the HPD interval must be <= this value
    # this gives us the tail of the distribution
    interval_idx_inc = int(np.floor(width * n))
    # number of values we are excluding
    # thus, possible number of HPD intervals
    # this gives us the head of the distribution
    n_intervals = n - interval_idx_inc
    # the width of each possible interval
    # for each possible head and tail value, 
    # what is the difference between them
    interval_width = [a_i - b_i for a_i, b_i in 
                      zip(dat[interval_idx_inc:], 
                          dat[:n_intervals])]
    # find the shortest interval
    min_idx = interval_width.index(min(interval_width))
    hpd_interval = (dat[min_idx], dat[min_idx+interval_idx_inc])
    dat_hpd = [item for item in dat if (item >= hpd_interval[0]) & (item <= hpd_interval[1])]
    dat_mid = np.quantile(dat_hpd, qs[1])
    return((hpd_interval[0], dat_mid, hpd_interval[1]))


def run():
    parser = argparse.ArgumentParser()
    # input files
    parser.add_argument('--metadata')
    parser.add_argument('--logFile')
    parser.add_argument('--plotParams', nargs=2, default=['seir.R0', 'time_of_mrca'])
    parser.add_argument('--absoluteTime', type=float)
    parser.add_argument('--timeFmt', default='day')
    parser.add_argument('--burnIn', default=0.10, type=float)
    parser.set_defaults(exclude=False)
    args = parser.parse_args()
    #args.logFile = 'data/SEIR_mu04_seed221110_unif10_32-52.log'
    #args.metadata = 't.tsv'

    #args.logFile = 'beast_comparison_wip/SEIR_simple_prop500_0928_seed123_100.log'
    if not args.absoluteTime:
        args.absoluteTime = pd.read_csv(args.metadata, sep='\t', header=None)[1].max()
    #args.absoluteTime = 0.26694
    #args.timeFmt = 'day'

    log = pd.read_csv(args.logFile, sep='\t', comment='#')
    # remove burnin
    log = log.iloc[int(log.shape[0]*args.burnIn):]

    # get time of mrca
    log['time_of_mrca'] = args.absoluteTime - log['TreeHeight']
    # convert if needed
    if args.timeFmt.lower()[:4] == 'day':
      # if year = 0, assume simualted data, assume 366 days in year like 2020
      if log['time_of_mrca'].apply(np.floor).min() <= 0:
        log['time_of_mrca'] *= 366
      else:
        log['time_of_mrca'] = log['time_of_mrca'].apply(datetime_from_numeric).apply(datetime.date.toordinal)

    print(log)
    # subset to just the plot data
    plot_dat = log[args.plotParams]
    plot_hpd = [hpd(plot_dat.iloc[:,0]), 
    	hpd(plot_dat.iloc[:,1])]


    # set up plot
    fig = plt.figure()
    gs = fig.add_gridspec(6, 6)
    ax1 = fig.add_subplot(gs[:2, :4])
    ax2 = fig.add_subplot(gs[2:, :4])
    ax3 = fig.add_subplot(gs[2:, 4:])

    # X axis density
    sns.kdeplot(data=plot_dat, x=args.plotParams[0], ax=ax1, color='#333333')
    # add vertical lines for 95% hpd
    ax1.axvline(plot_hpd[0][0], color='firebrick', ls='--')
    ax1.axvline(plot_hpd[0][2], color='firebrick', ls='--')
    # remove borders and ticks
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_xlabel(None)
    ax1.set_ylabel(None)
    [ax1.spines[i].set_visible(False) for i in ['bottom', 'top', 'left', 'right']]

    # Y axis density
    sns.kdeplot(data=plot_dat, y=args.plotParams[1], ax=ax3, color='#333333')
    # add vertical lines for 95% hpd
    ax3.axhline(plot_hpd[1][0], color='firebrick', ls='--')
    ax3.axhline(plot_hpd[1][2], color='firebrick', ls='--')
    # remove borders and ticks
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.set_xlabel(None)
    ax3.set_ylabel(None)
    [ax3.spines[i].set_visible(False) for i in ['bottom', 'top', 'left', 'right']]


    # scatter
    ax2.scatter(plot_dat.iloc[:,0], 
    	plot_dat.iloc[:,1],
    	alpha=0.1, 
    	color='#333333')

    # density
    sns.kdeplot(data=plot_dat, x=args.plotParams[0], y=args.plotParams[1], ax=ax2, 
    	levels=[0.05, 0.25, 0.5, 0.75, 0.95], color='firebrick', linestyles='--', bw_adjust=0.5)

    ax2.set_xlabel(r'$R_0$')
    ax2.set_ylabel('time of MRCA (days)')
    # need to make limits match
    ax1.set_xlim(ax2.get_xlim())
    ax3.set_ylim(ax2.get_ylim())


    fig.tight_layout()
    fig.savefig(f'{".".join(args.logFile.split(".")[:-1])}_{"_".join(args.plotParams)}.pdf')
    plt.close()



if __name__ == "__main__":
    run()
