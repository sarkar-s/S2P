#!/usr/bin/env python

# This file computes the average of pair correlation functions

import sys, os
from optparse import OptionParser
from multiprocessing import Pool
from itertools import repeat
import numpy as np
import math
import scipy.stats as spstats

def processOptions():
    parser = OptionParser()
    parser.add_option("-i", dest="datafolder", help="Name of the folder containing the outputs for a network growth simulation", default="")
    parser.add_option("-n", dest="samples", help="Number of network samples present in the network", default="")
    
    [options, args] = parser.parse_args()
    
    if options.datafolder != None:
        #os.chdir("..")
        os.chdir(str(options.datafolder))
    else:
        print 'Error wrong command line inputs'
        sys.exit()
        
    if options.samples!=None:
        return int(options.samples)
    else:
        print 'Error wrong command line inputs'
        sys.exit()

def extract_data_set(filename):
    all_lines = k_file.readlines()
    
    this_k, this_density = [], []
    
    for line in all_lines:
        this_set = line.rstrip('\r\n').split(',')
        
        this_k.append(float(this_set[0]))
        this_density.append(float(this_set[1]))
        
    return list(this_k), list(this_density)

if __name__ == '__main__':
    samples = processOptions()
    
    K_values, densities = [], []
    min_K, max_K = 1000.0, 0.0
    sample_size = 1000000
    
    good_file_count = 0
    
    for file_no in range(0,samples):
        try:
            k_file = open('k_distribution'+str(file_no)+'.csv','r')
            temp_K, temp_density = extract_data_set(k_file)
            k_file.close()
            
            K_values.append(temp_K)
            densities.append(temp_density)
            
            min_K = min(min_K,min(temp_K))
            max_K = max(max_K,max(temp_K))
            sample_size = min(sample_size,len(temp_K))

            good_file_count += 1
        except IOError:
            pass
        
    sample_size = int(sample_size)
    bin_size = (max_K - min_K)/float(sample_size-1)
    k_bins = np.linspace(min_K-0.5*bin_size,max_K+0.5*bin_size,sample_size+1)
    k_hist = np.zeros(sample_size)
    
    for kset in xrange(0,good_file_count):
        for i in xrange(0,len(K_values[kset])):
            for j in xrange(0,sample_size):
                if K_values[kset][i]>=k_bins[j] and K_values[kset][i]<k_bins[j+1]:
                    k_hist[j] += densities[kset][i]
                    break
    
    integral = sum(k_hist)*bin_size
    k_mids = []
    
    for i in xrange(0,sample_size):
        k_mids.append(0.5*(k_bins[i]+k_bins[i+1]))
        k_hist[i] *= 1.0/integral
        
    mean_k, var_k = 0.0, 0.0
    
    for i in xrange(0,sample_size):
        mean_k += k_mids[i]*k_hist[i]*bin_size
        var_k += (k_mids[i]**2)*k_hist[i]*bin_size
        
    var_k -= mean_k**2
    
    k_avg_file = open('averaged_k.csv','w')
    
    for i in xrange(0,len(k_mids)):
        print >> k_avg_file, k_mids[i],',',k_hist[i]
    
    k_avg_file.close()
    
    bin_size = 0.001
    
    x = np.linspace(bin_size/2,(1-bin_size/2),1000)
    
    # Compute gamma distribution parameters
    scale = var_k/mean_k
    shape = mean_k/scale
    
    print 'shape, scale =  ', shape, scale
    print 'mean, variance =  ', mean_k, var_k
    
    size = 1000
    
    x_b = np.linspace(spstats.gamma.ppf(0.0001,shape),spstats.gamma.ppf(0.9999,shape),size+1)
    
    x = []
    
    for i in xrange(0,size):
        x.append(0.5*(x_b[i]+x_b[i+1]))
    
    gamma_dist = spstats.gamma.pdf(x,shape)
    
    integral = sum(gamma_dist)*0.001
    
    #scale_factor = 1.0 #(max_K-min_K)/(x[-1]-x[0])
    
    for i in xrange(0,size):
        x[i] *= scale
    
    #x = x + min_K*np.ones(len(x))
    
    #smooth_k = sig.savgol_filter(k_hist,7,2)
    
    k_out_file = open('k_smooth.csv','w')
    
    #for i in xrange(0,len(k_mids)):
    #    print >> k_out_file, k_mids[i],',',smooth_k[i]
    
    gamma_dist *= 1.0/scale
    
    for i in xrange(0,len(x)):
        
        print >> k_out_file, x[i],',',gamma_dist[i]
    
    k_out_file.close()