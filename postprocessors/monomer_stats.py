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
    
    if options.datafolder != None:
        #os.chdir("..")
        os.chdir(str(options.datafolder))
    else:
        print 'Error wrong command line inputs'
        sys.exit()
        
    if options.samples!=None:
        return int(options.samples), float(options.time)
    else:
        print 'Error wrong command line inputs'
        sys.exit()

def extract_data_set(filename):
    data_lines = filename.readlines()
    
    # get headings from first line
    head_line = data_lines[0].rstrip('\r\n')
    labels = head_line.split(',')
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

def extract_chain_data(filename):
    data_lines = filename.readlines()
    
    # get headings from first line
    #head_line = data_lines[0].rstrip('\r\n')
    #labels = head_line.split(',')
    data_values = {}
    #data_x = {}
    
    #for name in labels:
    #    data_values[name] = []
    
    firstline = data_lines[1]
    lsize = len(firstline.rstrip('\r\n').split(','))
    
    data_values['Time'] = []
    data_values['DC'] = []
    labels = ['Time','DC']
    
    for i in xrange(1,lsize-1):
        name = 'Monomer'+str(i)
        data_values[name] = []
        labels.append(name)
    
    for line_no in range(1,len(data_lines)):
        data_set = data_lines[line_no].rstrip('\r\n').split(',')
        #data_x[labels[0]].append(float(data_set[0]))
        
        for k in data_values.keys():
            data_values[k].append(float(data_set[labels.index(k)]))
    
    return data_values

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

def get_intermediate_value(last_idx,this_time,time_set,data_set):
    data_value = data_set[last_idx]
    
    for idx in xrange(0,len(time_set)-1):
        if this_time>=time_set[idx] and this_time<=time_set[idx+1]:
            data_value = data_set[idx]
            last_idx = idx
            break
    
    return data_value, last_idx
                
if __name__ == '__main__':
    samples, calc_time = processOptions()
    good_file_count = 0
    
    timeDC_x = []
    timeDC_data = []
    sites_data = []
    
    max_time = 0.0
    min_time = 1000.0
    
    timeDC_mean = {}
    timeDC_max = {}
    timeDC_min = {}
    
    for file_no in range(0,samples):
        try:
            timedc_file = open('Time_DC'+str(file_no)+'.csv','r')
            
            monomerstats_file = open('monomer_count'+str(file_no)+'.csv','r')
            
            temp_x, temp_data = extract_data_set(timedc_file)
            max_time_for_this_set = max(temp_data['Time'])
            
            if max_time_for_this_set>=calc_time:
                chain_set = extract_chain_data(monomerstats_file)
                sites_data.append(chain_set)
                
                timeDC_x.append(temp_x)
                timeDC_data.append(temp_data)
                
                max_time = max(max_time,max(temp_data['Time']))
                min_time = min(min_time,min(temp_data['Time']))
                
                good_file_count += 1
            
            #timeDC_mean, timeDC_max, timeDC_min = process_data(timeDC_mean,timeDC_max,timeDC_min,timeDC_data,file_no)
            
            timedc_file.close()
            monomerstats_file.close()
        except IOError:
            pass
        
    time_limit = float(int(max_time))
    
    data_size = 1000
    
    l_t_set = np.linspace(np.log(min_time),np.log(time_limit),data_size)
    
    
    data_tags = []
    
    for key in sites_data[0].keys():
        if key not in ['Time','DC']:
            data_tags.append(key)
    
    #data_tags = ['Network Size','Molecule Diffusivity','Network Diffusivity']
        
    outfile = open('monomer_summary.csv','w')
    
    outstring = 'Time,DC'
    
    for key in data_tags:
        outstring += ','+key
    
    print >> outfile, outstring
    
    t_set = np.zeros(shape=(data_size+1))
    
    monomers = {}
    
    for key in data_tags:
        monomers[key] = np.zeros(shape=(data_size+1))
    #empty_mean = np.zeros(shape=(data_size+1))
    #surface_mean = np.zeros(shape=(data_size+1))
    #blocked_mean = np.zeros(shape=(data_size+1))
    dc_mean = np.zeros(shape=(data_size+1))

    t_idx = []
    
    for i in range(0,good_file_count):
        t_idx.append(int(0))

    for t_no in xrange(1,data_size+1):
        t_set[t_no] = math.exp(l_t_set[t_no-1])
        this_time = t_set[t_no]
        
        for k in range(0,good_file_count):
            dc_value, last_idx  = get_intermediate_value(t_idx[k],this_time,timeDC_data[k]['Time'],timeDC_data[k]['DC'])
            
            for tag in data_tags:
                this_value, last_idx = get_intermediate_value(t_idx[k],this_time,timeDC_data[k]['Time'],sites_data[k][tag])
                monomers[tag][t_no] += this_value
            
            #empty_value, last_idx = get_intermediate_value(t_idx[k],this_time,timeDC_data[k]['Time'],sites_data[k][data_tags[0]])
            #surface_value, last_idx = get_intermediate_value(t_idx[k],this_time,timeDC_data[k]['Time'],sites_data[k][data_tags[1]])
            #blocked_value, last_idx = get_intermediate_value(t_idx[k],this_time,timeDC_data[k]['Time'],sites_data[k][data_tags[2]])
            
            t_idx[k] = last_idx
            
            dc_mean[t_no] += dc_value
            
            #empty_mean[t_no] += empty_value
            #surface_mean[t_no] += surface_value
            #blocked_mean[t_no] += blocked_value
                
    for t_no in range(0,1001):
        outstring = str(t_set[t_no])+','+str(dc_mean[t_no]/good_file_count)
        
        for tag in data_tags:
            outstring += ','+str(monomers[tag][t_no]/good_file_count)
        
        print >> outfile, outstring
        
    outfile.close()  