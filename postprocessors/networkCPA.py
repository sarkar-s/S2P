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

class networkCPA(object):
    def __init__(self,_k_set,_k_hist,_k_mean):
        # No. of bins in the histogram
	self.K_mids = list(_k_set)
	self.K_hist = list(_k_hist)
        self.hist_size = len(_k_set)
	self.bin_size = self.K_mids[1] - self.K_mids[0]
	self.K_mean = _k_mean
        self.omega_max = 5.0*math.sqrt(self.K_mean)
        
	self.omega_set = list(np.linspace(0.0,self.omega_max,25))
        self.omega_set.pop(0)
        self.d_omega = self.omega_set[1] - self.omega_set[0]
        
        # Reduced Planck's constant
        #self.h_bar = (0.5/math.pi)*(6.626E-34)
        
        # Boltzmann factor
        #self.kB = 1.3805E-23
        #self.beta = 1.0/(self.kB*self.temperature)
                
        # K values range
        k_values = np.linspace(-math.pi,math.pi,51)
        
        # K integrator
        self.d_k = (k_values[1] - k_values[0])**3
        
        # K values range
        self.k_set = []
        
        for idx_k in range(0,8):
            self.k_set.append(0.5*(k_values[idx_k]+k_values[idx_k+1]))
        
        # T matrix epsilon
        self.epsilon = 1E-10
	
	# Tag for convergence
	self.convergence_tag = []

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

    def initialize_structs(self):
        """Initialize lists for the frequency dependent force constants
        and the frequencies at which they are evaluated using CPA.
        """
        self.Gamma = []
        self.z = []
        
        for omega in self.omega_set:
            self.z.append(complex(-(omega**2),(omega**2)*self.epsilon))
            self.Gamma.append(self.K_mean)
        
        # Brillouin Zone volume
        self.volume_BZ = 0.0
        
        for k_x in self.k_set:
            for k_y in self.k_set:
                for k_z in self.k_set:
                    self.volume_BZ += self.d_k

            
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
            print 'Gamma ', idx, ' calculated: ', self.Gamma[idx], self.count
            sys.stdout.flush()
	    if self.count==100:
		self.convergence_tag.append(False)
	    else:
		self.convergence_tag.append(True)
		
        #ofile = open('DoS_2.csv','w')
        
        #for idx in xrange(1,len(self.omega_set)):
        #    print >> ofile, self.omega_set[idx], ',', self.Gamma[idx]
        #
        #ofile.close()
        
    def get_this_gamma(self,idx):
        update = 1.0
        self.count = 0
        
        while (abs(update)/abs(self.Gamma[idx]))>1e-6 and self.count<100:
            update = self.integrate_with_distribution(idx) - self.Gamma[idx]
            self.Gamma[idx] += update
            self.count += 1
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

	#for idx in xrange(0,len(self.omega_set)):
	#    self.omega_set[idx] *= math.sqrt(self.mean_K)

    def write_DoS_MFP(self):
        self.phonon_file = open('phonon_values.csv','w')
        
        print >> self.phonon_file, 'Frequency',',','DOS',',','MFP',',','Wavelength',',','Velocity'
        
        for omega_idx in xrange(0,len(self.omega_set)):
	    if self.convergence_tag[omega_idx]==True:
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
        
        
        