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

class CoherentPA(object):
    def __init__(self,_hist_size,_omega_max,_temp):
        # No. of bins in the histogram
        self.hist_size = _hist_size
        self.omega_max = _omega_max
        self.temperature = _temp
        self.omega_set = list(np.linspace(0.0,self.omega_max,_hist_size))
        self.omega_set.pop(0)
        self.d_omega = self.omega_set[1] - self.omega_set[0]
        
        # Reduced Planck's constant
        self.h_bar = (0.5/math.pi)*(6.626E-34)
        
        # Boltzmann factor
        self.kB = 1.3805E-23
        self.beta = 1.0/(self.kB*self.temperature)
        
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
        self.local_stiffness_samples.append(local_stiffness)
        
        self.local_mass_samples.append(local_mass)

    def fitting_func(self,x,a,b,c):
	y = a + b*x + c*x**2
	
	return y

    def compute_histograms(self,true_global_mass,sampled_global_mass):
	mass_correction_factor = true_global_mass/float(sampled_global_mass)
	stiffness_samples = []
	
	for k in xrange(0,len(self.local_mass_samples)):
	    self.local_mass_samples[k] *= mass_correction_factor
	    
	    for j in xrange(0,6):
		stiffness_samples.append(self.local_stiffness_samples[k][j]/self.local_mass_samples[k])
		
	sample_size = len(stiffness_samples)
	self.mean_K = sum(stiffness_samples)/sample_size 

        max_stiffness = max(stiffness_samples)
	self.local_stiffness_samples = []
	
	for sample in stiffness_samples:
	    self.local_stiffness_samples.append(sample/self.mean_K)
        #self.local_stiffness_samples.append(0.0)

	# Compute interquartile range and histogram size
	#iqr = sp.stats.iqr(self.local_stiffness_samples)
	qrs = np.percentile(self.local_stiffness_samples,[75,25])
	self.min_K, self.max_K = min(self.local_stiffness_samples), max(self.local_stiffness_samples)
	hist_size = 0.5*(self.max_K-self.min_K)/(qrs[0]-qrs[1])
	hist_size *= len(self.local_stiffness_samples)**(1.0/3.0)
	self.hist_size = int(hist_size)
	
	#print 'Histogram calculation begins'
	#sys.stdout.flush()
	#
	#self.K_hist = {}
	#
	#for sample in self.local_stiffness_samples:
	#    this_key = int(1000*sample)/1000.0
	#    
	#    if this_key in self.K_hist.keys():
	#	self.K_hist[this_key] += 1.0/sample_size
	#    else:
	#	self.K_hist[this_key] = 1.0/sample_size

        # Histogram bins
        self.K_hist, self.K_bin_edges = np.histogram(self.local_stiffness_samples,self.hist_size,density=True)
        
        self.rho_hist, self.rho_bin_edges = np.histogram(self.local_mass_samples,self.hist_size/4,density=True)
        
        self.bin_size = self.K_bin_edges[1] - self.K_bin_edges[0]

	self.compute_mean_stiffness()

	#win_size = 3
	#
	#self.beta_dist = []
	#
	#for idx in xrange(0,self.hist_size):
	#    ll = max(idx - win_size,0)
	#    ul = min(idx+win_size,self.hist_size)
	#    
	#    pts_needed = (2*win_size+1) - (ul-ll)
	#
	#    if ll==0:
	#	ul += pts_needed
	#    else:
	#	ll -= pts_needed
	#
	#    popt, pcov = spopt.curve_fit(self.fitting_func,self.K_mids[ll:ul+1],self.K_hist[ll:ul+1])
	#
	#    self.beta_dist.append(self.fitting_func(self.K_mids[idx],popt[0],popt[1],popt[2]))
	
	# Create a smooth beta distribution of K_histogram
	#mu = np.mean(self.local_stiffness_samples)
	#sigma = np.std(self.local_stiffness_samples)
	#
	#beta_lambda = (mu-self.min_K)*(self.max_K-mu)/(sigma**2) - 1
	#
	#self.alpha = beta_lambda*(mu-self.min_K)/(self.max_K - self.min_K)
	#self.beta = beta_lambda*(self.max_K-mu)/(self.max_K - self.min_K)
	#
	#
	#
	#x = []
	#
	#for k in self.K_mids:
	#    x.append((k-self.min_K)/(self.max_K-self.min_K))
	#    
	#self.beta_dist = spstats.beta.pdf(np.array(x),self.alpha,self.beta)
	#
	#
	#print len(self.K_mids), len(self.beta_dist)
	#sys.stdout.flush()

    def write_histograms(self,proc_no):
        ofile_1 = open('K_hist'+str(proc_no)+'.csv','w')
        ofile_2 = open('rho_hist'+str(proc_no)+'.csv','w')
        
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
        
        # K_mean
#	for this_key in self.K_hist.keys():
#            self.K_mids.append(this_key)
#            self.K_mean += this_key*self.K_hist[this_key]

        for idx in xrange(0,self.hist_size-1):
            self.K_mids.append(0.5*(self.K_bin_edges[idx]+self.K_bin_edges[idx+1]))
            self.K_mean += self.K_mids[-1]*self.K_hist[idx]*(self.K_bin_edges[idx+1] - self.K_bin_edges[idx])
            
    def compute_modified_greens_function(self,z,omega,idx):
        self.tilde_G_z = 0.0
        
        for k_x in self.k_set:
            for k_y in self.k_set:
                for k_z in self.k_set:
                    integrand = math.cos(k_x) + math.cos(k_y) + math.cos(k_z)
                    self.tilde_G_z += self.d_k/((z/self.Gamma[idx]) + 6.0 - 2*integrand)
        
        self.tilde_G_z *= (1.0/(self.volume_BZ*self.Gamma[idx]))
        
        # Compute
        # [1 - z*tilde_G_z]/(3*Gamma_z)
        self.beta = (1.0 - z*self.tilde_G_z)/(3.0*self.Gamma[idx])
  
    def compute_Gamma_omega(self):
        for idx in xrange(0,len(self.omega_set)):
            self.get_this_gamma(idx)
            #print 'z: ', self.z[idx]
            #print 'Gamma ', idx, ' calculated: ', self.Gamma[idx]
            #sys.stdout.flush()
            
        #ofile = open('DoS_2.csv','w')
        
        #for idx in xrange(1,len(self.omega_set)):
        #    print >> ofile, self.omega_set[idx], ',', self.Gamma[idx]
        #
        #ofile.close()
        
    def get_this_gamma(self,idx):
        update = 1.0
        count = 0
        
        while (abs(update)/abs(self.Gamma[idx]))>0.001 and count<100:
            update = self.integrate_with_distribution(idx) - self.Gamma[idx]
            self.Gamma[idx] += update
            count += 1
            #print 'Update: ', abs(update)/abs(self.Gamma[idx]), self.Gamma[idx]
            #print self.tilde_G_z
            #sys.stdout.flush()
        
        #print count
    
    def integrate_with_distribution(self,idx):
        numerator = 0.0
        denominator = 0.0
        
        self.compute_modified_greens_function(self.z[idx],self.omega_set[idx],idx)
	
#	for this_key in self.K_hist.keys():
#            alpha = self.Gamma[idx] - this_key
#            factor = 1.0 - alpha*self.beta
#            
#            # Random variable
#            numerator += (this_key/factor)*self.K_hist[this_key]
#            denominator += (1.0/factor)*self.K_hist[this_key]

        for integ_idx in xrange(0,self.hist_size-1):
            alpha = self.Gamma[idx] - self.K_mids[integ_idx]
            factor = 1.0 - alpha*self.beta
            
            # Random variable
            numerator += (self.K_mids[integ_idx]/factor)*self.K_hist[integ_idx]*self.bin_size
            denominator += (1.0/factor)*self.K_hist[integ_idx]*self.bin_size
            
        update = numerator/denominator
            
        return update
    
    def check_integration(self):
        # Integration sum
        integ_sum = 0.0
        
        for integ_idx in xrange(0,self.hist_size-1):
            integ_sum += self.K_hist[integ_idx]*self.bin_size
        
        print 'Sum: ', integ_sum
        sys.stdout.flush()
        sys.exit()
    
    def compute_DoS_MFP(self):
        """This function computes the Density of States and Mean free path of the phonons
        from the complex frequency dependent stiffness constant.
        """
        self.g_omega = []
        self.mean_free_path = []
        self.sound_velocity = []
        self.wavelength = []
        
        for omega_idx in xrange(0,len(self.omega_set)):
            self.g_omega.append(2.0*math.pi*self.Gamma[omega_idx].imag/self.omega_set[omega_idx])
            self.sound_velocity.append(cm.sqrt(self.Gamma[omega_idx]))
            
            modulus = abs(self.Gamma[omega_idx])
            denominator = 2.0*self.omega_set[omega_idx]*self.sound_velocity[-1].imag
            
            self.mean_free_path.append(modulus/denominator)
            
            self.wavelength.append(2.0*math.pi*self.sound_velocity[-1].real/self.omega_set[omega_idx])

	for idx in xrange(0,len(self.omega_set)):
	    self.omega_set[idx] *= math.sqrt(self.mean_K)

    def write_DoS_MFP(self,proc_no):
        self.phonon_file = open('phonon_values'+str(proc_no)+'.csv','w')
        
        print >> self.phonon_file, 'Frequency',',','DOS',',','MFP',',','Wavelength',',','Velocity'
        
        for omega_idx in xrange(0,len(self.omega_set)):
            print >> self.phonon_file, self.omega_set[omega_idx],',',self.g_omega[omega_idx],',',self.mean_free_path[omega_idx],',',self.wavelength[omega_idx],',',self.sound_velocity[omega_idx].real
            
        self.phonon_file.close()
        
    def write_gamma_file(self,proc_no):
        self.gamma_file = open('gamma_values'+str(proc_no)+'.csv','w')
        
        for omega_idx in xrange(0,len(self.omega_set)):
            print >> self.gamma_file, self.omega_set[omega_idx],',',self.Gamma[omega_idx].real,',',self.Gamma[omega_idx].imag
            
        self.gamma_file.close()
        
    def compute_thermal_properties(self,proc_no):
        # Computing heat capacity and thermal conductivity
        self.heat_capacity = 0.0
        self.thermal_conductivity = 0.0
        
        for omega_idx in xrange(0,len(self.omega_set)):
            x = self.h_bar*self.true_omega[omega_idx]*self.beta
            
            print 'x: ', x
            sys.stdout.flush()
            
            integrand = (x**2)*math.exp(x)/((math.exp(x)-1)**2)
            self.heat_capacity += integrand*self.g_omega[omega_idx]
            self.thermal_conductivity += integrand*self.g_omega[omega_idx]*self.sound_velocity[omega_idx].real*self.mean_free_path[omega_idx]

        self.heat_capacity *= (self.kB*self.d_omega)
        self.thermal_conductivity *= (self.kB*self.d_omega)/3.0
        
        #ofile = open('phonon_properties'+str(proc_no)+'.csv','w')
        #
        #print >> ofile, self.heat_capacity
        #print >> ofile, self.thermal_conductivity
        #
        #ofile.close()
        
        
        