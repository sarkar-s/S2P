# Run gel modeling code
import sys, os
from optparse import OptionParser
from multiprocessing import Pool
from itertools import repeat
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
    head_line = data_lines[0].rstrip('\r\n')
    labels = head_line.split(',')
    data_values = {}
    
    for name in labels:
        data_values[name] = []
        
    for line_no in range(1,len(data_lines)):
        data_set = data_lines[line_no].rstrip('\r\n').split(',')
        for k in data_values.keys():
            data_values[k].append(float(data_set[labels.index(k)]))
    
    return labels[0], data_values

def extract_summary_data(filename):
    data_lines = filename.readlines()
    
    # Get summary keys
    summary_keys = data_lines[0].rstrip('\r\n').split(',')
    
    # Get summary values
    summary_value_strs = data_lines[1].rstrip('\r\n').split(',')
    
    summary_values = {}
    
    for idx in range(0,len(summary_value_strs)):
        summary_values[summary_keys[idx]] = float(summary_value_strs[idx])
    
    return summary_values

def process_summary_data(mean_set,max_set,min_set,data_set,file_no):
    if file_no==0:
        for key in data_set.keys():
            mean_set[key] = list(data_set[key])
            max_set[key] = list(data_set[key])
            min_set[key] = list(data_set[key])
    else:
        for key in data_set.keys():
            mean_set[key] = (1.0/float(file_no+1))*(file_no*mean_set[key] + data_set[key])
            max_set[key] = max(data_set[key],max_set[key])
            min_set[key] = min(data_set[key],min_set[key])
    
    return mean_set, max_set, min_set

def process_data(mean_set,max_set,min_set,var_set,data_set,file_no):
    if file_no==0:
        for key in data_set.keys():
            mean_set[key] = data_set[key]
            max_set[key] = data_set[key]
            min_set[key] = data_set[key]
    else:
        for key in data_set.keys():
            for id_no in range(0,len(data_set[key])):
                mean_set[key][id_no] += data_set[key][id_no]
                max_set[key][id_no] = max(data_set[key][id_no],max_set[key][id_no])
                min_set[key][id_no] = min(data_set[key][id_no],min_set[key][id_no])
    
    return mean_set, max_set, min_set
                
if __name__ == '__main__':
    samples = processOptions()
    
    all_means = {}
    all_var = {}
    
    all_values = {}
    
    for file_no in range(0,samples):
        summary_file = open('summary'+str(file_no)+'.csv','r')
        
        summary_data = extract_summary_data(summary_file)
        
        if file_no==0:
            for key in summary_data.keys():
                all_values[key] = [summary_data[key]]
        else:
            for key in summary_data.keys():
                all_values[key].append(summary_data[key])
           
        #summary_mean, summary_max, summary_min = process_summary_data(summary_mean,summary_max,summary_min,var_set,summary_data,file_no)
        
        summary_file.close()
        
    sum_file = open('total_summary.csv','w')
    
    for key in all_values.keys():
        all_means[key] = sum(all_values[key])/float(samples)
    
    for key in all_values.keys():
        variance = 0.0
        
        for samp_no in range(0,samples):
            variance += (all_values[key][samp_no] - all_means[key])**2
        
        variance *= (1.0/float(samples))
        
        all_var[key] = math.sqrt(variance)
    
    for key in summary_data.keys():
        print >> sum_file, key+','+str(all_means[key])+','+str(all_var[key])
    
    sum_file.close()