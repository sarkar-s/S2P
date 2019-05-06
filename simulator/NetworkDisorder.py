#!/usr/bin/env python

# Class for calculating density of states and mean free paths for the polymer network
import math
import sys
import numpy as np
import random as rand
import cmath as cm
import scipy.stats as spstats
#import scipy.signal as sig
import scipy.optimize as spopt

class NetworkDisorder(object):
    def __init__(self,):
        # Boltzmann factor
        self.kB = 1.3805E-23
        
        # List to contain samples of local stiffnesses and densities
        self.local_stiffness_samples = []
        self.local_mass_samples = []
        
        # K values range
        k_values = np.linspace(-math.pi,math.pi,51)
        
        # K integrator
        self.d_k = (k_values[1] - k_values[0])**3
        
        # K values range
        self.k_set = []
        
        for idx_k in range(0,50):
            self.k_set.append(0.5*(k_values[idx_k]+k_values[idx_k+1]))
        
        # T matrix epsilon
        self.epsilon = 1E-12

    def setup_scaling_factors(self,eps_self,lattice_h,molWt):
        # Compute frequency scaling factor
	wt_factor = molWt*(1.6605E-27)
	lattice_factor = (1E-20)*lattice_h**2
	energy_factor = (1.602E-19)*eps_self
	
	self.omega_factor = math.sqrt(energy_factor/(wt_factor*lattice_factor))

        factor_file = open('Scaling_factors.csv','w')
        
        print >> factor_file, 'Length factor',',',(lattice_h*(1E-10))
        print >> factor_file, 'Omega factor',',',self.omega_factor
        print >> factor_file, 'Mass factor',',',molWt*(1.6605E-27)

        factor_file.close()

    def add_stiffness_and_mass_to_bin(self,local_stiffness,local_mass):
        self.local_stiffness_samples.append(list(local_stiffness))
        
        self.local_mass_samples.append(local_mass)

    def compute_histograms(self,true_global_mass,sampled_global_mass):
	mass_correction_factor = true_global_mass/float(sampled_global_mass)
	stiffness_samples = []
	
	for k in xrange(0,len(self.local_mass_samples)):
	    self.local_mass_samples[k] *= mass_correction_factor
	    
	    for j in xrange(0,6):
		stiffness_samples.append(self.local_stiffness_samples[k][j])#/self.local_mass_samples[k])
		
	sample_size = len(stiffness_samples)
	self.mean_K = sum(stiffness_samples)/sample_size 

        max_stiffness = max(stiffness_samples)
        del self.local_stiffness_samples
	self.local_stiffness_samples = list(stiffness_samples)
	
	#for sample in stiffness_samples:
	#    self.local_stiffness_samples.append(sample)

	# Compute interquartile range and histogram size
	qrs = np.percentile(self.local_stiffness_samples,[75,25])
	self.min_K, self.max_K = min(self.local_stiffness_samples), max(self.local_stiffness_samples)
	hist_size = 0.5*(self.max_K-self.min_K)/(qrs[0]-qrs[1])
	hist_size *= len(self.local_stiffness_samples)**(1.0/3.0)
	self.K_hist_size = int(hist_size)

        # Histogram bins
        self.K_hist, self.K_bin_edges = np.histogram(self.local_stiffness_samples,self.K_hist_size,density=True)
        
        qrs = np.percentile(self.local_mass_samples,[75,25])
	self.min_m, self.max_m = min(self.local_mass_samples), max(self.local_mass_samples)
	hist_size = 0.5*(self.max_m-self.min_m)/(qrs[0]-qrs[1])
	hist_size *= len(self.local_mass_samples)**(1.0/3.0)
	self.m_hist_size = int(hist_size)

        self.rho_hist, self.rho_bin_edges = np.histogram(self.local_mass_samples,self.m_hist_size,density=True)
        
        self.bin_size = self.K_bin_edges[1] - self.K_bin_edges[0]

	self.compute_mean_stiffness()

    def write_histograms(self,proc_no):
        ofile_1 = open('K_distribution'+str(proc_no)+'.csv','w')
        ofile_2 = open('m_distribution'+str(proc_no)+'.csv','w')
        
#        for this_key in self.K_hist.keys():
#	    print >> ofile_1, this_key,',',self.K_hist[this_key]
	    
	for idx in xrange(0,len(self.K_mids)):
	    print >> ofile_1, self.K_mids[idx],',',self.K_hist[idx]

	for idx in xrange(0,len(self.rho_hist)):
            print >> ofile_2, 0.5*(self.rho_bin_edges[idx]+self.rho_bin_edges[idx+1]),',',self.rho_hist[idx]
        
        ofile_1.close()
        ofile_2.close()

    def initialize_structs(self):
        """Initialize lists for the frequency dependent force constants
        and the frequencies at which they are evaluated using CPA.
        """
        self.Gamma = []
        self.z = []
        
        for omega in self.omega_set:
            self.z.append(complex(-(omega**2),self.epsilon))
            self.Gamma.append(self.K_mean)
        
        # Brillouin Zone volume
        self.volume_BZ = 0.0
        
        for k_x in self.k_set:
            for k_y in self.k_set:
                for k_z in self.k_set:
                    self.volume_BZ += self.d_k
    
    def compute_mean_stiffness(self):
        """This function computes the average stiffness in the disordered network.
        """
        self.K_mids = []
        self.K_mean = 0.0

        for idx in xrange(0,self.K_hist_size-1):
            self.K_mids.append(0.5*(self.K_bin_edges[idx]+self.K_bin_edges[idx+1]))
            self.K_mean += self.K_mids[-1]*self.K_hist[idx]*(self.K_bin_edges[idx+1] - self.K_bin_edges[idx])
