#!/usr/bin/env python

# This file computes the average of pair correlation functions

import sys, os
from optparse import OptionParser
from multiprocessing import Pool
from itertools import repeat
import numpy as np
import math

def processOptions():
    parser = OptionParser()
    parser.add_option("-i", dest="datafolder", help="Name of the folder containing the outputs for a network growth simulation", default="")
    parser.add_option("-n", dest="samples", help="Number of network samples present in the network", default="")
    parser.add_option("-t", dest="time", help="Final simulation time", default="")
    
    [options, args] = parser.parse_args()
    
    if options.datafolder == None:
        print 'Error wrong command line inputs'
        sys.exit()
        
    if options.samples!=None:
        return str(options.datafolder), int(options.samples), float(int(options.time))
    else:
        print 'Error wrong command line inputs'
        sys.exit()
                
if __name__ == '__main__':
    datafolder, samples, calc_time = processOptions()
    
    current_dir = os.getcwd()
    
    # Compute Degree of Conversion statistics
    os.system('python DCstats.py -i '+str(datafolder)+' -n '+str(samples)+' -t '+str(calc_time))
    os.chdir(current_dir)
    
    # Compute Pair Correlation statistics
    os.system('python PCstats.py -i '+str(datafolder)+' -n '+str(samples))
    os.chdir(current_dir)
    
    #os.system('python PCvar_stats.py -i '+str(datafolder)+' -n '+str(samples))
    #os.chdir(current_dir)
    
    # Compute chain statistics
    os.system('python chainstats.py -i '+str(datafolder)+' -n '+str(samples)+' -t '+str(calc_time))
    os.chdir(current_dir)
    
    # Compute site statistics
    os.system('python sitestats.py -i '+str(datafolder)+' -n '+str(samples)+' -t '+str(calc_time))
    os.chdir(current_dir)
    
    # Compute diffusion statistics
    os.system('python diffusion_stats.py -i '+str(datafolder)+' -n '+str(samples)+' -t '+str(calc_time))
    os.chdir(current_dir)
    
    # Compute stiffness statistics
    os.system('python k_stats.py -i '+str(datafolder)+' -n '+str(samples))
    os.chdir(current_dir)
    
    # Compute stiffness statistics
    os.system('python monomer_stats.py -i '+str(datafolder)+' -n '+str(samples)+' -t '+str(calc_time))
    os.chdir(current_dir)
    
    # Compute stiffness statistics
    os.system('python movable_stats.py -i '+str(datafolder)+' -n '+str(samples)+' -t '+str(calc_time))
    os.chdir(current_dir)
    
    # Compute fluctuation statistics
    #os.system('python fluctuationstats.py -i '+str(datafolder)+' -n '+str(samples)+' -t '+str(calc_time))
    #os.chdir(current_dir)
    
    # Compute phonon statistics
    #os.system('python phonon_stats.py -i '+str(datafolder)+' -n '+str(samples))
    #os.chdir(current_dir)
    
    # Compute density statistics
    #os.system('python rho_stats.py -i '+str(datafolder)+' -n '+str(samples))
    #os.chdir(current_dir)
    
    # Compute summary statistics
    #os.system('python SummaryStats.py -i '+str(datafolder)+' -n '+str(samples))
    #os.chdir(current_dir)