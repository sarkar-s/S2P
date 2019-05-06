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
    data_lines = filename.readlines()
    
    # get headings from first line
    #head_line = data_lines[0].rstrip('\r\n')
    #labels = head_line.split(',')
    labels = ['rho','prob mass']
    data_values = {}
    data_x = {}
    
    for name in labels:
        if name==labels[0]:
            data_x[name] = []
        else:
            data_values[name] = []
        
    for line_no in range(0,len(data_lines)):
        data_set = data_lines[line_no].rstrip('\r\n').split(',')
        data_x[labels[0]].append(float(data_set[0]))
        
        for k in data_values.keys():
            data_values[k].append(float(data_set[labels.index(k)]))
    
    return data_x, data_values

def process_data(mean_set,max_set,min_set,data_set,file_no):
    if file_no==0:
        for key in data_set.keys():
            mean_set[key] = list(data_set[key])
            max_set[key] = list(data_set[key])
            min_set[key] = list(data_set[key])
    else:
        for key in data_set.keys():
            for id_no in range(0,len(data_set[key])):
                x = mean_set[key][id_no]
                mean_set[key][id_no] = (1.0/float(file_no+1))*(float(file_no)*x + data_set[key][id_no])
                y = max_set[key][id_no]
                max_set[key][id_no] = max(data_set[key][id_no],y)
                z = min_set[key][id_no]
                min_set[key][id_no] = min(data_set[key][id_no],z)
    
    return mean_set, max_set, min_set

def get_rho_value(low_rho,up_rho,this_set,mass_set):
    prob_mass = 0.0
    
    for idx in xrange(0,len(this_set)):
        if this_set[idx]>=low_rho and this_set[idx]<up_rho:
            prob_mass += mass_set[idx]

    return prob_mass

if __name__ == '__main__':
    samples = processOptions()
    good_file_count = 0
    
    rho_x = []
    rho_data = []
    
    max_rho = 0.0
    min_rho = 1000.0
    
    for file_no in range(0,samples):
        try:
            rho_file = open('rho_hist'+str(file_no)+'.csv','r')
            
            temp_x, temp_data = extract_data_set(rho_file)
            
            rho_x.append(temp_x)
            rho_data.append(temp_data)
            
            max_rho = max(max_rho,max(temp_x['rho']))
            min_rho = min(min_rho,min(temp_x['rho']))
            
            good_file_count += 1
            
            rho_file.close()
        except IOError:
            pass
        
    data_size = 40

    rho_set = np.linspace(min_rho,max_rho,data_size)
    rho_set_length = len(rho_set)
    
    rho_mids = []
    
    for rho_no in xrange(0,rho_set_length-1):
        rho_mids.append(0.5*(rho_set[rho_no]+rho_set[rho_no]))
    
    outfile = open('rho_summary.csv','w')
        
    for rho_no in range(0,rho_set_length-1):
        mass_mean = 0.0

        for k in range(0,good_file_count):
            mass_value = get_rho_value(rho_set[rho_no],rho_set[rho_no+1],rho_x[k]['rho'],rho_data[k]['prob mass'])
            mass_mean += mass_value
        
        print >> outfile, str(rho_mids[rho_no])+','+str(mass_mean/good_file_count)
        
    outfile.close()  