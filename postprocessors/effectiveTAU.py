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
import matplotlib.pyplot as plt
from matplotlib import rc

class effectiveTAU(object):
    def __init__(self,_k_set,_k_hist,_k_mean,_bins,_omega_max,_omega_min):
        # No. of bins in the histogram
	self.K_mids = list(_k_set)
	self.K_hist = list(_k_hist)
        self.hist_size = len(_k_set)
	self.bin_size = self.K_mids[1] - self.K_mids[0]
	self.bins =_bins
	self.K_mean = _k_mean
        self.omega_max = _omega_max
	self.omega_min = _omega_min

	self.omega_set = list(np.linspace(self.omega_min,self.omega_max,251))
        #self.omega_set.pop(0)
        self.d_omega = self.omega_set[1] - self.omega_set[0]
        
        # Reduced Planck's constant
        #self.h_bar = (0.5/math.pi)*(6.626E-34)
        
        # Boltzmann factor
        #self.kB = 1.3805E-23
        #self.beta = 1.0/(self.kB*self.temperature)
                
        # K values range
        k_values = np.linspace(0.0,math.pi,11)
        
        # K integrator
        self.d_k = (k_values[1] - k_values[0])**3
        
        # K values range
        self.k_set = []
        
        for idx_k in range(0,10):
            self.k_set.append(0.5*(k_values[idx_k]+k_values[idx_k+1]))
        
        # T matrix epsilon
        self.epsilon = 1e-0

	self.k_xi = np.linspace(0.005,0.995,100)
	self.dk = self.k_xi[1]-self.k_xi[0]

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
	self.tilde_G_z = []
	self.tilde_lambda = []
        self.convergence = []

	guess = self.compute_initial_guess()

        for omega in self.omega_set:
	    epsilon = self.K_mean*(omega**2)/self.epsilon
	    
            self.z.append(complex(-(omega**2),self.epsilon))
            self.Gamma.append(complex(guess,0.0))
	    self.tilde_G_z.append(0.0)
	    self.tilde_lambda.append(0.0)
            self.convergence.append(False)

        # Brillouin Zone volume
        self.volume_BZ = 0.0
        
        for k_x in self.k_set:
            for k_y in self.k_set:
                for k_z in self.k_set:
                    self.volume_BZ += self.d_k

    def compute_initial_guess(self):
	initial_Gamma = self.K_mean
	update = initial_Gamma
	new_gamma = 0.0
	
	while abs(update/initial_Gamma)>(1e-6):
	    new_gamma = 0.0
	    for i in range(0,self.hist_size):
		new_gamma += self.bin_size*self.K_mids[i]/(1+(1.0/3.0)*(self.K_mids[i]/initial_Gamma - 1))
		
	    update = new_gamma - initial_Gamma
	    initial_Gamma = new_gamma
	    
	print self.K_mean, initial_Gamma
	    
	return initial_Gamma

    def compute_modified_greens_function(self,z,omega,idx):
        self.tilde_G_z[idx] = 0.0
        
        for k_x in self.k_set:
            for k_y in self.k_set:
                for k_z in self.k_set:
                    integrand = math.cos(k_x) + math.cos(k_y) + math.cos(k_z)
                    self.tilde_G_z[idx] += self.d_k/(z + self.Gamma[idx]*(6.0 - 2*integrand))
        
        self.tilde_G_z[idx] *= (1.0/(self.volume_BZ))
	
	#for k in self.k_xi:
	#    self.tilde_G_z[idx] += 3*(k**2)*(1.0/(z + (k**2)*self.Gamma[idx]))*self.dk
	
        # Compute
        # [1 - z*tilde_G_z]/(3*Gamma_z)
        self.beta = (1.0 - z*self.tilde_G_z[idx])/(3.0*self.Gamma[idx])
    
    def compute_lambda(self,z,idx):
        self.tilde_lambda[idx] = (1.0 - z*self.tilde_G_z[idx])/self.Gamma[idx]
  
    def compute_Gamma_omega(self):
	idset = range(0,len(self.omega_set))
	#idset.reverse()
        #for idx in xrange(0,len(self.omega_set)):
	for idx in idset:
            self.get_this_gamma(idx)
	    if self.count<1000:
		self.convergence[idx] = True
            #print 'z: ', self.z[idx]
            print 'Gamma ', idx, ' calculated: ', self.z[idx]*self.tilde_G_z[idx], self.z[idx], self.count
            sys.stdout.flush()
	    if idx<idset[-1]:
		self.Gamma[idx+1] = self.Gamma[idx]

        #ofile = open('DoS_2.csv','w')
        
        #for idx in xrange(1,len(self.omega_set)):
        #    print >> ofile, self.omega_set[idx], ',', self.Gamma[idx]
        #
        #ofile.close()
        
    def get_this_gamma(self,idx):
        update = complex(1.0,1.0)
        self.count = 0
        
        #while (abs(update.imag)/abs(self.Gamma[idx].imag))>1e-9 and self.count<1000:
	while (abs(update.imag))>1e-12 and self.count<1000:
            update = self.integrate_with_distribution(idx)# - self.Gamma[idx]
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
	self.compute_lambda(self.z[idx],idx)
	
#	for this_key in self.K_hist.keys():
#            alpha = self.Gamma[idx] - this_key
#            factor = 1.0 - alpha*self.beta
#            
#            # Random variable
#            numerator += (this_key/factor)*self.K_hist[this_key]
#            denominator += (1.0/factor)*self.K_hist[this_key]

        for integ_idx in xrange(0,self.hist_size):
            #alpha = self.Gamma[idx] - self.K_mids[integ_idx]
            #factor = 1.0 - alpha*self.beta
	    factor = 1.0 + (1.0/3.0)*(self.K_mids[integ_idx]-self.Gamma[idx])*self.tilde_lambda[idx]
	#    if self.omega_set[idx]<math.sqrt(self.K_mean):
	#	factor = 1.0 + (1.0/3.0)*(self.K_mids[integ_idx]-self.Gamma[idx])/self.Gamma[idx]
	#    else:
	#	factor = 1.0 + (1.0/3.0)*(self.K_mids[integ_idx]-self.Gamma[idx])*self.tilde_lambda[idx]
	    
            # Random variable
            numerator += ((self.K_mids[integ_idx]-self.Gamma[idx])/factor)*self.K_hist[integ_idx]*self.bins[integ_idx]
            denominator += (1.0/factor)*self.K_hist[integ_idx]*self.bins[integ_idx]
            
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

	for idx in range(len(self.omega_set)-1,-1,-1):
	    if self.convergence[idx]==False:
                self.convergence.pop(idx)
                self.omega_set.pop(idx)
                self.Gamma.pop(idx)
                self.tilde_G_z.pop(idx)
        
        for omega_idx in xrange(0,len(self.omega_set)):
            self.g_omega.append(-2.0*self.tilde_G_z[omega_idx].imag/(math.pi*self.omega_set[omega_idx]))
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
    
    def plot_results(self):
	"""This function plots the results.
	"""
	plt.plot(self.omega_set,self.g_omega,marker='o',ms=7.5,lw=0.0)
	plt.minorticks_on()
	plt.grid()
	plt.xlim(0.0,max(self.omega_set))
	plt.xticks(size=20)
	plt.yticks(size=20)
	plt.xlabel(r'Frequency, $(\omega)$',size=20)
	plt.ylabel(r'Reduced DOS, $g(\omega)/\omega^2$',size=20)
	plt.show()

        
        
        