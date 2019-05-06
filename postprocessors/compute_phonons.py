#!/usr/bin/env python

import sys, os
from optparse import OptionParser
#from networkCPA import *
import numpy as np
import scipy.stats as spstats
import scipy.signal as sig
from latticeCPA import *
from continuumCPA import *

def processOptions():
    parser = OptionParser()
    parser.add_option("-i", dest="datafolder", help="Name of the folder containing the outputs for a network growth simulation", default="")
    #parser.add_option("-n", dest="samples", help="Number of network samples present in the network", default="")
    
    [options, args] = parser.parse_args()
    
    if options.datafolder != None:
        #os.chdir("..")
        os.chdir(str(options.datafolder))
    else:
        print 'Error wrong command line inputs'
        sys.exit()
        
    #if options.samples!=None:
    #    return int(options.samples)
    #else:
    #    print 'Error wrong command line inputs'
    #    sys.exit()
        
def extract_data_set(k_file):
    all_lines = k_file.readlines()
    
    this_k, this_density = [], []
    
    for line in all_lines:
        this_set = line.rstrip('\r\n').split(',')
        
        this_k.append(float(this_set[0]))
        this_density.append(float(this_set[1]))
        
    return list(this_k), list(this_density)
        
if __name__ == '__main__':
    processOptions()
    #
    #K_values, densities = [], []
    #min_K, max_K = 1000.0, 0.0
    #sample_size = 1000000
    #
    #for file_no in range(0,samples):
    #    k_file = open('K_distribution'+str(file_no)+'.csv','r')
    #    temp_K, temp_density = extract_data_set(k_file)
    #    k_file.close()
    #    
    #    K_values.append(temp_K)
    #    densities.append(temp_density)
    #    
    #    min_K = min(min_K,min(temp_K))
    #    max_K = max(max_K,max(temp_K))
    #    sample_size = min(sample_size,len(temp_K))
    #
    #sample_size = int(sample_size)
    #bin_size = (max_K - min_K)/float(sample_size-1)
    #
    #k_bins = np.linspace(min_K-0.5*bin_size,max_K+0.5*bin_size,sample_size+1)
    #
    #k_hist = np.zeros(sample_size)
    #
    ##for bin_no in xrange(0,sample_size):
    #for kset in xrange(0,samples):
    #    for i in xrange(0,len(K_values[kset])):
    #        for j in xrange(0,sample_size):
    #            if K_values[kset][i]>=k_bins[j] and K_values[kset][i]<k_bins[j+1]:
    #                k_hist[j] += densities[kset][i]
    #                break
    #
    #integral = sum(k_hist)*bin_size
    #k_mids = []
    #
    #for i in xrange(0,sample_size):
    #    k_mids.append(0.5*(k_bins[i]+k_bins[i+1]))
    #    k_hist[i] *= 1.0/integral
    #    
    #mean_k, var_k = 0.0, 0.0
    #
    #for i in xrange(0,sample_size):
    #    mean_k += k_mids[i]*k_hist[i]*bin_size
    #    var_k += (k_mids[i]**2)*k_hist[i]*bin_size
    #    
    #var_k -= mean_k**2
    
    ksmooth_file = open('k_smooth.csv','r')
    klines = ksmooth_file.readlines()
    x, dist, bins = [], [], []
    mean_x = 0.0
    var_x = 0.0
    
    for line in klines:
        this_set = line.rstrip('\r\n').split(',')
        x.append(float(this_set[0]))
        dist.append(float(this_set[1]))
        mean_x += x[-1]*dist[-1]
        var_x += (x[-1]**2)*dist[-1]
        
        if len(x)>1:
            bins.append(x[-1]-x[-2])
    
    mean_x *= bins[-1]
    var_x *= bins[-1]
    var_x += -mean_x*mean_x
    
    bins.append(bins[-1])
    
    _omega_max = 5*math.sqrt(mean_x)
    _omega_min = 0.1
    
    scale = var_x/mean_x
    shape = mean_x/scale
    
    print mean_x, var_x
    print 'shape, scale = ',shape, scale
    
    #phonon_calculator = continuumCPA(x,dist,mean_x,bins,1.5,0.1)
    phonon_calculator = latticeCPA(x,dist,mean_x,bins,_omega_max,_omega_min)
    phonon_calculator.initialize_structs()
    phonon_calculator.compute_Gamma_omega()
    phonon_calculator.compute_DoS_MFP()
    phonon_calculator.write_DoS_MFP()
    phonon_calculator.plot_results()