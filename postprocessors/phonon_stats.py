#!/usr/bin/env python

# This file computes the average of phonon calculations

import sys, os
from optparse import OptionParser
from multiprocessing import Pool
from itertools import repeat
import numpy as np

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
    head_line = data_lines[0].rstrip('\r\n')
    labels = head_line.split(',')
    for idx in range(0,len(labels)):
        labels[idx] = labels[idx].rstrip(' ').lstrip(' ')
        
    print labels
    data_values = {}
    data_x = {}
    
    for name in labels:
        if name==labels[0]:
            data_x[name] = []
        else:
            data_values[name] = []
        
    for line_no in range(1,len(data_lines)):
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

def get_density_mass(w_low,w_up,w_set,mass_set,path_set,wave_set,vel_set):
    total_mass = 0.0
    total_path = 0.0
    total_wave = 0.0
    total_vel = 0.0
    
    for w_no in xrange(0,len(w_set)):
        if w_set[w_no]>=w_low and w_set[w_no]<w_up:
            total_mass += mass_set[w_no]
            total_path += mass_set[w_no]*path_set[w_no]
            total_wave += mass_set[w_no]*wave_set[w_no]
            total_vel += mass_set[w_no]*vel_set[w_no]

    return total_mass, total_path, total_wave, total_vel

if __name__ == '__main__':
    samples = processOptions()
    good_files_count = 0
    
    pc_mean = {}
    pc_max = {}
    pc_min = {}
    all_pc = []
    
    for file_no in range(0,samples):
        try:
            correlation_file = open('phonon_values'+str(file_no)+'.csv','r')
            
            pc_x, pc_data = extract_data_set(correlation_file)
            pc_mean, pc_max, pc_min = process_data(pc_mean,pc_max,pc_min,pc_data,file_no)
            
            all_pc.append(pc_data)
            
            good_files_count += 1
            
            correlation_file.close()
        except IOError:
            pass
        
    outfile = open('phonon_values_summary.csv','w')
    
    data_size = 80
    
    w_set = np.linspace(0.0,max(pc_x['Frequency']),data_size)
    
    w_mids = []
    
    for i in xrange(0,data_size-1):
        w_mids.append(0.5*(w_set[i]+w_set[i+1]))
    
    x_key = pc_x.keys()[0]
    
    for w_no in xrange(0,data_size-1):
        mass_mean = 0.0
        path_mean = 0.0
        wave_mean = 0.0
        vel_mean = 0.0

        for k in range(0,good_files_count):
            mass_value, path_value, wave_value, vel_value = get_density_mass(w_set[w_no],w_set[w_no+1],pc_x['Frequency'],all_pc[k]['DOS'],all_pc[k]['MFP'],all_pc[k]['Wavelength'],all_pc[k]['Velocity'])
            mass_mean += mass_value
            path_mean += path_value
            vel_mean += vel_value
            wave_mean += wave_value
        
        path_mean *= 1.0/mass_mean
        wave_mean *= 1.0/mass_mean
        vel_mean *= 1.0/mass_mean
    
        print >> outfile, str(w_mids[w_no])+','+str(mass_mean/good_files_count)+','+str(path_mean/good_files_count)+','+str(wave_mean/good_files_count)+','+str(vel_mean/good_files_count)
        
        
    outfile.close()  