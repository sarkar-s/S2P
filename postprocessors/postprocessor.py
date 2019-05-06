# Run gel modeling code
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
        os.chdir("./"+options.datafolder)
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
            mean_set[key] = data_set[key]
            max_set[key] = data_set[key]
            min_set[key] = data_set[key]
    else:
        for key in data_set.keys():
            mean_set[key] = (1.0/float(file_no+1))*(file_no*mean_set[key] + data_set[key])
            max_set[key] = max(data_set[key],max_set[key])
            min_set[key] = min(data_set[key],min_set[key])
    
    return mean_set, max_set, min_set

def process_data(mean_set,max_set,min_set,data_set,file_no):
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
    
    timeDC_mean, pc_mean, summary_mean = {}, {}, {}
    timeDC_max, pc_max, summary_max = {}, {}, {}
    timeDC_min, pc_min, summary_min = {}, {}, {}
    
    for file_no in range(0,samples):
        #timeDC_file = open('Time_DC'+str(file_no)+'.csv','r')
        #correlation_file = open('pair_correlations'+str(file_no)+'.csv','r')
        summary_file = open('summary'+str(file_no)+'.csv','r')
        
        #timeDC_x, timeDC_data = extract_data_set(timeDC_file)
        #timeDC_mean, timeDC_max, timeDC_min = process_data(timeDC_mean,timeDC_max,timeDC_min,timeDC_data,file_no)
        #
        #pc_x, pc_data = extract_data_set(correlation_file)
        #pc_mean, pc_max, pc_min = process_data(timeDC_mean,timeDC_max,timeDC_min,timeDC_data,file_no)
        
        summary_data = extract_summary_data(summary_file)
        summary_mean, summary_max, summary_min = process_summary_data(summary_mean,summary_max,summary_min,summary_data,file_no)
        
        #timeDC_file.close()
        #correlation_file.close()
        summary_file.close()
        
    for key in summary_data.keys():
        print summary_mean[key], summary_max[key], summary_min[key]
    
        