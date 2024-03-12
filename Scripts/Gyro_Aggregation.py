#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 09:08:04 2020

@author: davidclark
"""


from __future__ import division
import pandas as pd
import os
import numpy as np
import multiprocessing as mp
import scipy.stats as stats
import macroeco.models as mod
import traceback

# Check for pypartitions
try:
    import pypartitions as pyp # Package obtained from https://github.com/klocey/partitions
except ImportError:
    print("Feasible set package not found. Download from\n" +
                "https://github.com/klocey/partitions and set to PYTHONPATH")

"""
Description
------------

Script for generating feasible set predictions of log variance for the data 
described in Johnson and Wilber, using the algorithms provided in the 
pypartitions package from Locey and McGlinn 2013.  

The work horse function is `feasible_mixture` which wraps around
pypartitions `rand_partitions`. One can just use either of these functions to
generate random partitions for a given P and H.

"""

def multiprocess_feasible_tpl(processes, ph_data, samples, max_P=10000):
    """
    A function for multiprocessing to speed up the calculations
    """

    pool = mp.Pool(processes=processes)

    results = [pool.apply_async(calc_tpl_var,
            args=(row_id, P, H, samples, max_P)) for row_id, P, H in ph_data]

    results = [p.get() for p in results]
    results.sort() # to sort the results
    return results

def calc_tpl_var(row_id, P, H, samples, max_P=10000):
    """
    Draw feasible set. Use MaxEnt approximation if P > 10000. Return summary 
    statistics from feasible set draw.

    Description
    -----------
    row_id : int
        Integer used for logging purposed.  Does not affect analysis
    P : int
        Total number of parasites
    H : int
        Total number of hosts
    samples : int
        Number of macrostates to draw form the feasible set
    max_P : int
        Above this value a MaxEnt approximation is used to speed up calculations.
    """

    print("Working on sample {0} with, P={1}, H={2}".format(row_id, P, H))
    if P < max_P:
        # Use feasible set
        vals = feasible_mixture([(P, H)], samples=samples, center="median")[0]
    else:
        #pred_var = mod.cnbinom(mu=P / H, k_agg=1, b=P).var()
        vals = mod.cnbinom(mu=P / H, k_agg=1, b=P).rvs(samples*H).reshape(samples, H)

    # Extract both untransformed and log10-transformed statistics
    pred_var = np.var(vals, ddof=1, axis=1)
    uncert = np.var(pred_var, ddof=1)
    log_uncert = np.var(np.log10(pred_var), ddof=1)
    lower, median, upper = stats.scoreatpercentile(pred_var, (2.5, 50.0, 97.5))
    log_lower, log_median, log_upper = stats.scoreatpercentile(np.log10(pred_var), (2.5, 50.0, 97.5))
    mean = np.mean(pred_var)
    log_mean = np.mean(np.log10(pred_var))

    return((row_id, median, uncert, lower, upper, mean, 
                    log_median, log_uncert, log_lower, log_upper, log_mean)) 

def feasible_mixture(para_host_vec, samples=200, center="mean"):
    """
    Gives a feasible set mixture from a given parasite host vector.

    Parameters
    -----------
    para_host_vec : list of tuples
        Length of the list is the number of heterogeneities. Each tuple is
        (P, H) where P is parasites and H is hosts. e.g. [(100, 10), (20, 30)]

    samples : int
        Number of feasible set samples to take

    center : str
        Either "mean", "median", or "mode".  They give very similar answers.  The mean
        will guarantee that the mean of the returned predicted vector will
        equal the expected mean for all sample sizes.  This is not necessarily
        True for the median or mode for small sample sizes. For the mode, if
        multiple values are found for the mode the minimum value is taken. The
        mode measure is more dependent on sample size than the mean or median,
        though it is what is used in Locey and White 2013.

    Returns
    -------
    : 2D array, 1D array
        The sorted array of all the sampled feasible sets
        The predicted center of the feasible set ("mean", "median", "mode")

    Examples
    --------
    full_feas, cent_feas = feasible_mixture([(100, 10), (20, 30)],
                                samples=200,  center="median")

    """

    mixed_sample = []
    para_host_vec = convert_to_int(para_host_vec)

    for ph in para_host_vec:

        if ph[0] == 0: # Hosts don't have any parasites

            tfeas = np.zeros(ph[1] * samples).reshape(samples, ph[1])
            mixed_sample.append(tfeas)

        else:

            tfeas = pyp.rand_partitions(int(ph[0]), int(ph[1]), samples, zeros=True)
            mixed_sample.append(tfeas)

    mix_feas = np.concatenate(mixed_sample, axis=1)
    sorted_feas = np.sort(mix_feas)[:, ::-1]

    if center == "mean":
        med_feas = np.mean(sorted_feas, axis=0)
    elif center == "median":
        med_feas = np.median(sorted_feas, axis=0)
    else:
        med_feas = stats.mode(sorted_feas, axis=0)[0].flatten()

    return sorted_feas, med_feas


def convert_to_int(para_host_vec):
    """
    Takes in a list of tuples and makes sure every item is a integer
    """

    para, host = zip(*para_host_vec)
    H = np.sum(host)
    para_round = np.round(para, decimals=0)
    host_round = np.round(host, decimals=0)

    # Make sure the vector still adds up. Should only miss by one (need to test)
    if np.sum(host) < H:
        ind = np.argmax(np.array(host) - np.floor(host))
        hosts_round[ind] = host_round[ind] + 1

    elif np.sum(host) > H:
        ind = np.argmin(np.array(host) - np.floor(host))
        hosts_round[ind] = host_round[ind] - 1


    return(zip(para_round.astype(np.int), host_round.astype(np.int)))

if __name__ == '__main__':

    # Load TPL data
    dat = pd.read_csv("HP_Aggregation_MJ_males.csv")

    tdat = dat.sort_values(by="P")
    tdat.index.name = "index"
    ph_data = zip(tdat.index, tdat.P, tdat.H)

    # Calculate feasible set variances using 1000 samples and join them to the table
    res = multiprocess_feasible_tpl(3, ph_data, 1000)
    res_df = pd.DataFrame(res)
    res_df.rename(columns={0: "index", 1:"fs_var", 2:"uncert", 3:"lower", 
                           4:"upper", 5:"fs_var_mean", 6:"fs_var_log_med", 
                           7:"uncert_log", 8:"lower_log", 
                           9:"upper_log", 10:"fs_var_mean_log"}, inplace=True)
    res_df.set_index("index", inplace=True)
    jdat = tdat.join(res_df)

    # Save the results of the analysis.  Will be the same as the previous 
    # file with slightly different estimates due to sampling error.
    jdat.to_csv("tpl_data_MJ_males.csv", index=False)


