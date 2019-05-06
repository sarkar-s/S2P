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
    head_line = data_lines[0].rstrip('\r\n')
    labels = head_line.split(',')
    data_values = {}
    #data_x = {}
    
    for name in labels:
        data_values[name] = []
        
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
    chain_data = []
    
    max_time = 0.0
    min_time = 1000.0
    
    timeDC_mean = {}
    timeDC_max = {}
    timeDC_min = {}
    
    for file_no in range(0,samples):
        try:
            timedc_file = open('Time_DC'+str(file_no)+'.csv','r')
            
            chainstats_file = open('chainstat'+str(file_no)+'.csv','r')
            
            temp_x, temp_data = extract_data_set(timedc_file)
            max_time_for_this_set = max(temp_data['Time'])
            
            if max_time_for_this_set>=calc_time:
                chain_set = extract_chain_data(chainstats_file)
                chain_data.append(chain_set)
                
                timeDC_x.append(temp_x)
                timeDC_data.append(temp_data)
                
                max_time = max(max_time,max(temp_data['Time']))
                min_time = min(min_time,min(temp_data['Time']))
                
                good_file_count += 1
            
            #timeDC_mean, timeDC_max, timeDC_min = process_data(timeDC_mean,timeDC_max,timeDC_min,timeDC_data,file_no)
            
            timedc_file.close()
            chainstats_file.close()
        except IOError:
            pass
        
    time_limit = float(int(max_time))
    
    data_size = 1000
    
    l_t_set = np.linspace(np.log(min_time),np.log(time_limit),data_size)
    
    data_tag = 'Active chains'
        
    outfile = open('active_summary.csv','w')
    
    t_set = np.zeros(shape=(data_size+1))
    jn_mean = np.zeros(shape=(data_size+1))
    jn_max = np.zeros(shape=(data_size+1))
    jn_min = 100.0*np.ones(shape=(data_size+1))

    t_idx = []
    
    for i in range(0,good_file_count):
        t_idx.append(int(0))

    for t_no in xrange(1,data_size+1):
        t_set[t_no] = math.exp(l_t_set[t_no-1])
        this_time = t_set[t_no]
        
        for k in range(0,good_file_count):
            data_value, last_idx = get_intermediate_value(t_idx[k],this_time,timeDC_data[k]['Time'],chain_data[k][data_tag])
            
            t_idx[k] = last_idx
            
            jn_mean[t_no] += data_value
            jn_max[t_no] = max(data_value,jn_max[t_no])
            jn_min[t_no] = min(data_value,jn_min[t_no])
            
            # If this time is greater than equal to
            #if t_idx[k]==0 and this_time<(timeDC_data[k]['Time'][int(t_idx[k]+1)]):
            #    jn_mean[t_no] += 0.0
            #    jn_max[t_no] = 0.0
            #    jn_min[t_no] = 0.0
            #elif timeDC_data[k]['Time'][int(t_idx[k])]==timeDC_data[k]['Time'][-1]:
            #    jn_mean[t_no] += chain_data[k][data_tag][int(t_idx[k])]
            #    jn_max[t_no] = max(jn_max[t_no],chain_data[k][data_tag][int(t_idx[k])])
            #    jn_min[t_no] = min(jn_min[t_no],chain_data[k][data_tag][int(t_idx[k])])
            #elif this_time>=timeDC_data[k]['Time'][int(t_idx[k])] and this_time<(timeDC_data[k]['Time'][int(t_idx[k]+1)]):
            #    jn_mean[t_no] += chain_data[k][data_tag][int(t_idx[k])]
            #    jn_max[t_no] = max(jn_max[t_no],chain_data[k][data_tag][int(t_idx[k])])
            #    jn_min[t_no] = min(jn_min[t_no],chain_data[k][data_tag][int(t_idx[k])])
            #elif this_time>=timeDC_data[k]['Time'][int(t_idx[k])] and this_time>=(timeDC_data[k]['Time'][int(t_idx[k]+1)]):
            #    t_idx[k] += 1
            #    jn_mean[t_no] += chain_data[k][data_tag][int(t_idx[k])]
            #    jn_max[t_no] = max(jn_max[t_no],chain_data[k][data_tag][int(t_idx[k])])
            #    jn_min[t_no] = min(jn_min[t_no],chain_data[k][data_tag][int(t_idx[k])])
    
    for t_no in range(0,1001):
        print >> outfile, str(t_set[t_no])+','+str(jn_mean[t_no]/good_file_count)+','+str(jn_max[t_no]-jn_min[t_no])
        
    outfile.close()  