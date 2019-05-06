#!/usr/bin/env python

import sys, os
from optparse import OptionParser
#from networkCPA import *
import numpy as np
import scipy.stats as spstats
import scipy.signal as sig
from phonon_properties import *
from continuumCPA import *
import math

def processOptions():
    parser = OptionParser()
    parser.add_option("-i", dest="datafolder", help="Name of the folder containing the outputs for a network growth simulation", default="")
    parser.add_option("-s", dest="scale", help="Scaling factor for the frequency axis", default="1")
    
    [options, args] = parser.parse_args()
    
    if options.datafolder != None:
        #os.chdir("..")
        os.chdir(str(options.datafolder))
        
        return float(options.scale)
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

def compute_DB_value():
    max_x = 1.0
    x = np.linspace(0.0,1.0/max_x,1001)
    x_mids = []
    
    for i in xrange(0,1000):
        x_mids.append(0.5*(x[i]+x[i+1]))

    integral = 0.0
    
    for i in xrange(0,1000): 
        factor = (math.exp(x_mids[i])/(math.exp(x_mids[i])-1)**2)*x_mids[i]**4
        
        integral += factor*(x[i+1]-x[i])
        
    integral *= 3*max_x**3
        
    print integral
  
if __name__ == '__main__':
    scale = processOptions()
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
    
    #compute_DB_value()
    
    ksmooth_file = open('phonon_values.csv','r')
    klines = ksmooth_file.readlines()
    omega, g_omega, MFP, lambda_omega, v_omega = [], [], [], [], []
    ksmooth_file.close()
    
    mean_x = 0.0
    var_x = 0.0
    
    for line in klines[1:]:
        this_set = line.rstrip('\r\n').split(',')
        
        omega.append(float(this_set[0]))
        g_omega.append(float(this_set[1]))
        MFP.append(float(this_set[2]))
        lambda_omega.append(float(this_set[3]))
        v_omega.append(float(this_set[4]))
    
    phonon_properties = phonon_properties(omega,g_omega,MFP,lambda_omega,v_omega,scale)
    phonon_properties.compute_cv_kappa()
    phonon_properties.write_results()