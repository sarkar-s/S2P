import os, sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from optparse import OptionParser

def processOptions():
    parser = OptionParser()
    parser.add_option("-f", dest="inputfile", help="Name of the file containing the inputs for a photon transport simulation.", default="")
    
    [options, args] = parser.parse_args()
    
    return options.inputfile
    

input_file = processOptions()

rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

data_lines = open(input_file,'r').readlines()

y_label = input_file.rstrip('.txt')

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

x_key = labels[0]
x_arr = np.array(data_values[labels[0]])

for i in range(1,len(labels)):
    lval = 1.0*i
    plt.plot(x_arr,np.array(data_values[labels[i]]),lw=lval,label=labels[i])

plt.minorticks_on()
plt.grid()
plt.xticks(size='x-large')
plt.yticks(size='x-large')
plt.xlabel(labels[0],size='x-large')
plt.ylabel(y_label,size='x-large')
plt.legend(loc='best',fontsize='x-large')
plt.show()

