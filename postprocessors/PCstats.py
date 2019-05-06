#!/usr/bin/env python

# This file computes the average of pair correlation functions

import sys, os
from optparse import OptionParser
from multiprocessing import Pool
from itertools import repeat

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
                
if __name__ == '__main__':
    samples = processOptions()
    good_file_count = 0
    
    pc_mean = {}
    
    for file_no in range(0,samples):
        try:
            correlation_file = open('pc_final_'+str(file_no)+'.csv','r')
            pc_x, pc_data = extract_data_set(correlation_file)
            
            if pc_mean=={}:
                for key in pc_data.keys():
                    pc_mean[key] = []
                    pc_mean[key].append(pc_data[key])
            else:
                for key in pc_data.keys():
                    pc_mean[key].append(pc_data[key])
            
            #pc_mean.append(pc_data)
            
            correlation_file.close()
            good_file_count += 1
        except IOError:
            pass
        
    outfile = open('pc_final_summary.csv','w')
    
    x_key = pc_x.keys()[0]
    
    for id_no in range(0,len(pc_x[x_key])):
        outstring = str(pc_x[x_key][id_no])
        for key in pc_data.keys():
            total_value = 0
            for k in range(0,good_file_count):
                total_value += pc_mean[key][k][id_no]
                
            outstring += ','+str(total_value/good_file_count)
            
        print >> outfile, outstring
        
    outfile.close()
    
    good_file_count = 0
    
    pc_mean = {}
    
    for file_no in range(0,samples):
        try:
            correlation_file = open('pc_light_off_'+str(file_no)+'.csv','r')
            
            pc_x, pc_data = extract_data_set(correlation_file)
            
            if pc_mean=={}:
                for key in pc_data.keys():
                    pc_mean[key] = []
                    pc_mean[key].append(pc_data[key])
            else:
                for key in pc_data.keys():
                    pc_mean[key].append(pc_data[key])
            
            correlation_file.close()
            good_file_count += 1
        except IOError:
            good_file_count += -1
            pass
        
    outfile = open('pc_light_off_summary.csv','w')
    
    x_key = pc_x.keys()[0]
    
    for id_no in range(0,len(pc_x[x_key])):
        outstring = str(pc_x[x_key][id_no])
        for key in pc_data.keys():
            total_value = 0
            for k in range(0,good_file_count):
                total_value += pc_mean[key][k][id_no]
                
            outstring += ','+str(total_value/good_file_count)
            
        print >> outfile, outstring
        
    outfile.close()  