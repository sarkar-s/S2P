# Gel kinetics model for cross linking polymers

import os, sys
import numpy as np
import math
import cmath

class EffectiveDiffusion(object):
    def __init__(self,input_dict):
#        # Interaction energy between the monomer units
#        self.eps_self = float(input_dict['Self interaction energy with initiators (in eV)'])
#
#	# Interaction energy between initiators and monomers
#	self.eps_I = float(input_dict['Interaction energy with initiators (in eV)'])
#	
#	# Interaction energy between co-initiators and monomers
#	self.eps_co_I = float(input_dict['Interaction energy with initiators (in eV)'])
#	
#	# Interaction energy with other monomer molecules
#        self.s_eps = float(input_dict['Bond angle stiffness (in eV)'])
        # Lattice size
        self.lattice_size = float(input_dict['Lattice half size'])
	self.temp = float(input_dict['Temperature (in K)'])

	# Boltzmann factor(kBT) in eV
	self.beta = (1.0/0.0259)*(300.0/self.temp)

	# Total sites
	self.total_sites = 8*self.lattice_size**3
	
	# Total node sites
	self.total_node_sites = self.lattice_size**3
	
	self.phi = np.zeros(3)
	self.phi[0] = 1.0

    def set_energy_levels(self,eps_quanta):
	self.eps = np.zeros(6)
	
	self.eps[0] = 0.0
	self.eps[1] = eps_quanta/2.0
	self.eps[2] = eps_quanta
	self.eps[3] = eps_quanta
	self.eps[4] = 1.5*eps_quanta
	self.eps[5] = 2.0*eps_quanta
	
	#self.eps[0] = 0.0
	#self.eps[1] = eps_quanta
	#self.eps[2] = eps_quanta
	#self.eps[3] = 2.0*eps_quanta
	#self.eps[4] = 2.0*eps_quanta
	#self.eps[5] = 2.0*eps_quanta
	
	#self.eps[0] = 0.0
	#self.eps[1] = eps_quanta/float(8.0)
	#self.eps[2] = eps_quanta/float(4.0)
	#self.eps[3] = 3.0*eps_quanta/float(8.0)
	#self.eps[4] = eps_quanta/float(2.0)
	#self.eps[5] = 5.0*eps_quanta/float(8.0)
	#self.eps[6] = 3.0*eps_quanta/float(4.0)
	#self.eps[7] = eps_quanta
	
	self.g = []
	
	for i in xrange(0,6):
	    self.g.append(math.exp(-self.beta*self.eps[i]))

    def compute_interface_probabilities(self,j_sites,j_neighbor_sites,p_sites,p_neighbor_sites):
	self.phi = np.zeros(3)
	
	self.phi[2] = (j_neighbor_sites)/float(self.total_sites)
	self.phi[1] = (p_neighbor_sites - j_neighbor_sites)/float(self.total_sites)
	self.phi[0] = 1.0 - sum(self.phi[1:])
	
	#self.phi[4] = j_sites/float(self.total_sites)
	#self.phi[3] = j_neighbor_sites/float(self.total_sites)
	#self.phi[2] = (p_sites - j_sites)/float(self.total_sites)
	#self.phi[1] = (p_neighbor_sites - j_neighbor_sites)/float(self.total_sites)
	#self.phi[0] = 1.0 - sum(self.phi[1:])
	
	self.pop_out_string = str(self.phi[0])
	
	for i in xrange(1,3):
	    self.pop_out_string += ','+str(self.phi[i])
	    
	self.pop_out_string += ','+str(p_neighbor_sites/float(self.total_sites))+','+str(j_neighbor_sites/float(self.total_sites))
	
	if j_sites>0:
	    self.pop_out_string += ','+str(j_neighbor_sites/float(j_sites))
	
    def compute_mean_conductivities(self):
	self.g_bar = np.zeros(3)
	
	self.g_bar[0] = self.phi[0]*self.g[0] + self.phi[1]*self.g[1] + self.phi[2]*self.g[3]
	self.g_bar[1] = self.phi[0]*self.g[1] + self.phi[1]*self.g[2] + self.phi[2]*self.g[4]
	self.g_bar[2] = self.phi[0]*self.g[3] + self.phi[1]*self.g[4] + self.phi[2]*self.g[5]
	
	#if (self.phi[0]+self.phi[1]+self.phi[3])>0.0:
	#    self.g_bar[0] = (self.phi[0]+self.g[1]*self.phi[1] + self.g[4]*self.phi[3])/(self.phi[0]+self.phi[1]+self.phi[3])
	#else:
	#    self.g_bar[0] = 0
	#
	#self.g_bar[1] = (self.g[1]*self.phi[0] + self.g[2]*self.phi[1] + self.g[3]*self.phi[2] + self.g[5]*self.phi[3])
	#    
	#if (self.phi[1]+self.phi[3])>0.0:
	#    self.g_bar[2] = (self.g[3]*self.phi[1] + self.g[6]*self.phi[3])/(self.phi[1]+self.phi[3])
	#else:
	#    self.g_bar[2] = 0.0
	#    
	#self.g_bar[3] = (self.g[4]*self.phi[0] + self.g[5]*self.phi[1] + self.g[6]*self.phi[2] + self.g[7]*self.phi[3])
	
	self.g_bar_string = str(self.g_bar[0])
	
	for i in xrange(1,3):
	    self.g_bar_string += ','+str(self.g_bar[i])
	    
    def set_initial_guess(self,_guess):
	self.ini_g_em = _guess

    def compute_effective_diffusivity(self,this_eps):
	self.set_energy_levels(this_eps)
	self.compute_mean_conductivities()

	#g_em = np.inner(self.g_bar,self.phi)
	
	g_em = self.ini_g_em
	
	count = 0
	correction = 1.0
	
	while abs(correction)>1e-3:
	    count += 1
	    #fn = self.compute_f_g_em(g_em)
	    #fprime = self.compute_f_prime_g_em(g_em)
	    #
	    #delta_g_em = -fn/fprime
	    
	    numer, denom = self.compute_coeffs(g_em)
	    
	    delta_g_em = -numer/denom
	
	    g_em += delta_g_em
	    
	    correction = delta_g_em/g_em
	
	return g_em

    def compute_coeffs(self,g_em):
	numer, denom = 0.0, 0.0
	
	for i in xrange(0,3):
	    numer += self.phi[i]*(g_em - self.g_bar[i])/(2*g_em + self.g_bar[i])
	    denom += self.phi[i]/(2*g_em + self.g_bar[i])
	    
	return numer, denom

    def compute_f_g_em(self,g_em):
	fn = 0.0
	
	for i in xrange(0,3):
	    fn += self.phi[i]*(g_em - self.g_bar[i])/(2*g_em + self.g_bar[i])
	
	return fn

    def compute_f_prime_g_em(self,g_em):
	fprime = 0.0
	
	for i in xrange(0,3):
	    fprime += self.phi[i]/(2*g_em + self.g_bar[i])
	    fprime += -2*self.phi[i]*(g_em - self.g_bar[i])/((2*g_em + self.g_bar[i])**2)
	    
	return fprime