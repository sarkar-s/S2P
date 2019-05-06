# Gel kinetics model for cross linking polymers

import os, sys
import numpy as np
import math
import cmath

class EffectiveDiffusion_camboni(object):
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

	# Boltzmann factor(kBT) in eV
	self.beta = (1.0/0.0259)

	# Total sites
	self.total_sites = 8*self.lattice_size**3
	
	# Total node sites
	self.total_node_sites = self.lattice_size**3
    
    def set_energy_levels(self,eps_quanta):
	self.eps_crosslinks = 6.0*eps_quanta
	self.eps_node = 3.0*eps_quanta

    def compute_interface_probabilities(self,blocked_sites,neighbor_sites):
	self.update_population(blocked_sites,neighbor_sites)

    def compute_camboni_diffusivity(self,this_eps):
	if self.phi3<=(2.0/3.0):
	    self.compute_coeffs(this_eps)
	    D_em = self.compute_roots()
	else:
	    D_em = 0.0
	
	return D_em

    def compute_coeffs(self,this_eps):
	self.E_mean = self.phi1 + self.phi2*math.exp(-0.5*self.beta*this_eps)
	self.E_bar = (1.0/3.0)*(1.0 + math.exp(-0.5*self.beta*this_eps))
	
	self.a = 1.0
	self.b = self.E_mean*(3.0*self.E_bar - 3.0*self.E_mean)
	self.c = (self.E_mean**2)*(1.0 - 3.0*(1.0 - self.phi3))*math.exp(-0.5*self.beta*this_eps)
	
	self.root_factor = 1.0/((3.0-1.0)*(self.phi1 + self.phi2*math.exp(-self.beta*this_eps)))

    def compute_roots(self):
	disc = math.sqrt(self.b**2 - 4.0*self.a*self.c)
	
	self.r1 = self.root_factor*(-self.b+disc)/(2.0*self.a)
	self.r2 = self.root_factor*(-self.b-disc)/(2.0*self.a*self.root_factor)
	
	if self.r1>=0 and self.r1<=1.0:
	    f_em = self.r1
	else:
	    f_em = self.r2
	    
	return f_em

    def compute_effective_diffusivity(self,this_eps):#,junction_sites):
	#self.update_population(blocked_sites,neighbor_sites,surface_sites)#,junction_sites)
	
	self.compute_camboni_conductivities(this_eps)
		
	if self.phi3<=(2.0/3.0):
	    self.compute_camboni_coeffs()
	    f_em = self.compute_camboni_roots()
	else:
	    f_em = 0.0
	    
	return f_em

    def compute_camboni_conductivities(self,this_eps):
	self.f1 = 0.0
	self.f2 = 0.0
	
	self.f1 = self.phi1 + self.phi2*math.exp(-0.5*self.beta*this_eps)
	self.f2 = self.phi1*math.exp(-0.5*self.beta*this_eps)+self.phi2*math.exp(-self.beta*this_eps)
    
    def compute_camboni_coeffs(self):
	self.a = 4.0
	self.b = -self.f1*(6.0*self.phi1-2.0)-self.f2*(6.0*self.phi2-2.0)
	self.c = -self.f1*self.f2*(2.0-3.0*self.phi3)
	
    def compute_camboni_roots(self):
	disc = math.sqrt(self.b**2 - 4.0*self.a*self.c)
	
	self.r1 = (-self.b+disc)/(2.0*self.a)
	self.r2 = (-self.b-disc)/(2.0*self.a)
	
	if self.r1>=0 and self.r1<=1.0:
	    f_em = self.r1
	else:
	    f_em = self.r2
	    
	return f_em
	
    # Computes coefficient to the polynomial
    #def compute_coeffs(self):
    #    self.E_avg = self.phi1 + self.phi2*math.exp(-self.beta_i*self.eps)+self.phi3*math.exp(-0.5*self.beta_i*self.eps)
    #    
    #    self.c4 = 4.0
    #    self.c3 = 6.0*self.E_bar - 4.0*self.E_avg - 2.0*self.phi1
    #    self.c2 = math.exp(-1.5*self.beta_i*self.eps)*(1.0 - 3.0*(self.phi2+self.phi3))
    #    self.c2 += math.exp(-self.beta_i*self.eps)*(1.0 - 3.0*(self.phi1+self.phi2))
    #    self.c2 += math.exp(-0.5*self.beta_i*self.eps)*(1.0 - 3.0*(self.phi1+self.phi3))
    #    self.c1 = -math.exp(-1.5*self.beta_i*self.eps)*(1.0-1.5*self.phi4)

    def function(self,f_trial):
        func_val = self.c4*(f_trial**3)
        func_val += self.c3*(f_trial**2)
        func_val += self.c2*f_trial
        func_val += self.c1

    def f_prime(self,f_trial):
	func_prime = 3*self.c4*(f_trial**2)
	func_prime += 2*self.c3*f_trial
	func_prime += self.c2
        
    # Diffusivity function
    def get_roots(self):
	d0 = self.c3**2 - 3.0*self.c4*self.c2
	d1 = 2.0*self.c3**3 - 9.0*self.c4*self.c3*self.c2 + 27.0*(self.c4**2)*self.c1

	C_t = ((d1 + cmath.sqrt(d1**2 - 4.0*d0**3))/2.0)**(1.0/3.0)
	factor1 = -(1.0/(3.0*self.c4))
	self.x1 = factor1*(self.c3 + C_t + d0/C_t)
	self.x2 = factor1*(self.c3 + self.r2*C_t + d0/(self.r2*C_t))
	self.x3 = factor1*(self.c3 + self.r3*C_t + d0/(self.r3*C_t))

	#print 'Roots : ', self.x1, self.x2, self.x3

    def get_root(self):
	self.get_roots()
	
	rx1 = self.x1.real
	rx2 = self.x2.real
	rx3 = self.x3.real

	if rx1 > 0.0:
		self.f_em = rx1
	elif rx2 > 0.0:
		self.f_em = rx2
	else:
		self.f_em = rx3

	return self.f_em

    def get_diffusion_factor(self,empty_sites,diffusible_interfaces):
	D_factor = (1.0/48.0)*diffusible_interfaces/empty_sites
	
	return D_factor
    
    def update_population(self,blocked_sites,neighbor_sites):
	total_blocked_sites = blocked_sites
	
	self.diffusible_sites = neighbor_sites
	
	self.empty_sites = self.total_sites - (blocked_sites + neighbor_sites)
	
	self.phi3 = total_blocked_sites/self.total_sites
	self.phi2 = neighbor_sites/self.total_sites
	self.phi1 = 1.0 - self.phi2 - self.phi3
	
	#self.phi4 = (neighbor_sites - surface_sites)/self.total_sites
	
	self.empty_sites *= (1.0/self.total_sites)
