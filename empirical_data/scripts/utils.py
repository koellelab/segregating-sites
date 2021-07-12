# move to utils
def format_seq_name(seq, sep, field, remove):
    description = seq.description
    for i in remove:
        description = description.replace(i, '')
    description = description.split(sep)[field]
    seq.description = ''
    seq.name = ''
    seq.id = description
    return(seq)


def convert_to_datetime(dt, fmt):
    # pandas to_datetime infers the day and month when it's missing
    # we want to return NAN when that is the case
    from datetime import datetime
    import numpy as np
    try: 
        return(datetime.strptime(dt, fmt))
    except:
        return(np.nan)

        
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
    days_elapsed = int(((numdate%1)+1e-10)*days_in_year) -1
    date = datetime.datetime(int(numdate),1,1) + datetime.timedelta(days=days_elapsed)
    return date


# move to utils
# from treetime
# https://github.com/neherlab/treetime/blob/1e378bfbbb98451ce4d6dd309d3699323102ee64/treetime/utils.py
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


# move to utils
def matlab_ordinal_from_numeric(numdate):
    import datetime
    # takes numeric date (e.g. 2020.01) as input
    # and returns MATLAB ordinal date
    year = int(numdate)
    ordinal_year = \
        datetime.datetime.strptime(f'{year}', '%Y').toordinal()
    from calendar import isleap
    days_in_year = 366 if isleap(year) else 365
    # adjust by 1 because ordinal year includes the first day of the year
    days_elapsed = numdate%1 * days_in_year - 1
    # 366 adjustment for matlab
    ordinal_date = ordinal_year + days_elapsed + 366
    return(ordinal_date)




