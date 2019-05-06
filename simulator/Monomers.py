#!/usr/bin/env python

import numpy as np
import random as rand
import math
import sys

class Monomers(object):
    def __init__(self,input_dict):
        """This class stores the following information about the monomers in the matrix
        Mole fraction
        Functionality
        """
        self.no_monomers = int(input_dict['Monomer types'])
        self.occupancy = float(input_dict['Site occupancy fraction'])
        
        # Mole fraction of monomers
        self.monomer_fraction = []
        
        # Functionality of monomers
        self.functionality = []
        
        # Molecular weight of monomers
        self.molecular_wt = []
        
        # Molecular volume of monomers
        self.molecular_vol = []
        
        # Molecular hopping rates
        self.hopping_rate = []
        
        # Interaction energy of monomers
        self.eps_self = []
        self.mean_eps_self = 0.0
        
        # Diffusivity factors
        self.d_factor_monomers = []
        
        # Propagation rate constant of monomers
        self.k_p = np.zeros((self.no_monomers,self.no_monomers))
        
        # Termination rate constant of monomers
        self.k_t = np.zeros((self.no_monomers,self.no_monomers))
        
        # Termination rate constant of monomers
        self.eps = np.zeros((self.no_monomers,self.no_monomers))
                
        # Probability distribution of monomers
        self.prob_set = np.zeros(shape=(self.no_monomers))
        
        # Population of monomers added to the network
        self.added_set = np.zeros(shape=(self.no_monomers))
        
        # Network size
        self.network_size = 0
        
        for mo_type in range(0,self.no_monomers):
            prefix = 'Monomer '+str(mo_type + 1)+': '
            
            self.monomer_fraction.append(float(input_dict[prefix+'Mole fraction']))
            self.functionality.append(float(input_dict[prefix+'Functionality']))
            self.molecular_wt.append(float(input_dict[prefix+'Molecular weight (in Da)']))
            self.molecular_vol.append(float(input_dict[prefix+'Molecular volume (in Angstroms^3)']))
            self.eps_self.append(float(input_dict[prefix+'Self interaction energy (in eV)']))
            self.hopping_rate.append(float(input_dict[prefix+'Hopping rate']))
            self.d_factor_monomers.append(1.0)
            
        self.mean_eps_self = sum(self.eps_self)/float(self.no_monomers)
        
        rate_constants_file = input_dict['Rate constants file']
        self.read_rate_constants(rate_constants_file)
        
        term_constants_file = input_dict['Termination constants file']
        self.read_term_constants(rate_constants_file)
        
        #interactions_file = input_dict['Interaction energy file']
        #self.read_interactions(interactions_file)
        
        self.interaction_range = int(input_dict['Interaction range'])
        
        self.create_interaction_matrix()
    
    def read_rate_constants(self,rate_constants_file):
        file_handle = open(rate_constants_file, 'r')
        rate_lines = file_handle.readlines()
        file_handle.close()

        rr = 0
        
        for line in rate_lines:
            cc = 0
            rate_set = line.rstrip('\r\n').split('\t')
            
            for value in rate_set:
                self.k_p[rr,cc] = float(value)
                cc += 1
            
            rr += 1
    
    def read_term_constants(self,term_constants_file):
        file_handle = open(term_constants_file, 'r')
        term_lines = file_handle.readlines()
        file_handle.close()
        
        rr = 0
        
        for line in term_lines:
            cc = 0
            term_set = line.rstrip('\r\n').split('\t')
            
            for value in term_set:
                self.k_t[rr,cc] = float(value)
                cc += 1
            
            rr += 1
            
    def read_interactions(self,interactions_file):
        file_handle = open(interactions_file, 'r')
        term_lines = file_handle.readlines()
        file_handle.close()
        
        rr = 0
        
        for line in term_lines:
            cc = 0
            term_set = line.rstrip('\r\n').split('\t')
            
            for value in term_set:
                self.eps[rr,cc] = float(value)
                cc += 1
            
            rr += 1
            
    def create_interaction_matrix(self):
        self.intermat = np.zeros(shape=(self.no_monomers,self.no_monomers))
        
        for i in xrange(0,self.no_monomers):
            for j in xrange(0,self.no_monomers):
                self.intermat[i,j] = math.sqrt(self.eps_self[i]*self.eps_self[j])
    
    def compute_monomer_populations(self,density,lattice_half_size,initiator_sites):
        """This class computes the total monomer population based on the lattice size
        """
        self.total_monomers = self.total_effective_units - initiator_sites
        self.monomer_population = []
        self.initial_monomers = []
        self.total_sites = lattice_half_size**3
        
        self.total_functionality = 0.0
        
        for mo_type in range(0,self.no_monomers):
            self.monomer_population.append(int(self.monomer_fraction[mo_type]*self.total_monomers))
            self.initial_monomers.append(self.monomer_population[-1])
            self.total_functionality += self.monomer_population[-1]*self.functionality[mo_type]
        
        self.monomer_population.append(self.total_monomers - sum(self.monomer_population))
        
        self.free_monomer_fraction = self.total_monomers/float(self.total_sites)
        
        self.construct_prob_dist()
        
    def update_monomer_population(self,monomer_type,update_type='remove'):
        """Updates the free monomers available after each polymerization event
        """
        if update_type=='remove':
            self.monomer_population[monomer_type] += -1
            self.added_set[monomer_type] += 1
            
            self.total_monomers += -1
            self.network_size += 1
        elif update_type=='add':
            self.monomer_population[monomer_type] += 1
            self.added_set[monomer_type] += -1
            
            self.total_monomers += 1
        
        self.free_monomer_fraction = self.total_monomers/float(self.total_sites)

        self.construct_prob_dist()
        self.compute_mean_eps_self()
        
    def construct_prob_dist(self):
        total_prob_mass = 0.0
        
        for mo_type in range(0,self.no_monomers):
            self.prob_set[mo_type] = float(self.monomer_population[mo_type])*self.d_factor_monomers[mo_type]
            total_prob_mass += self.prob_set[mo_type]
            
        self.mean_d_factors = total_prob_mass/float(self.total_monomers)
        
        for mo_type in range(0,self.no_monomers):
            self.prob_set[mo_type] *= 1.0/total_prob_mass
    
    def compute_mean_eps_self(self):
        self.mean_eps_self = 0.0
        
        for mo_type_1 in range(0,self.no_monomers):
            total_eps = 0.0
        
            for mo_type_2 in range(0,self.no_monomers):
                effective_eps = math.sqrt(self.eps_self[mo_type_1]*self.eps_self[mo_type_2])
                total_eps += effective_eps*float(self.initial_monomers[mo_type_2]-self.monomer_population[mo_type_2])
            
            self.mean_eps_self += total_eps/float(self.network_size)
        
        self.mean_eps_self *= 1.0/float(self.no_monomers)
    
    def select_a_monomer(self):
        xi = rand.random()
        total = 0.0
        
        #selected_monomer = 0
        
        for mo_type in range(0,self.no_monomers):
            if total<=xi and xi<(total+self.prob_set[mo_type]):
                selected_monomer = mo_type
            
            total += self.prob_set[mo_type]
        
        return selected_monomer
    
    def compute_domain_volume(self,density,lattice_size_half):
        #mass_factor = 1.66054E-24
        #
        self.total_sites = lattice_size_half**3
        #
        self.mean_molecular_vol = 0.0
        self.mean_molecular_wt = 0.0
        
        for mo_type in xrange(0,self.no_monomers):
            self.mean_molecular_vol += self.monomer_fraction[mo_type]*self.molecular_vol[mo_type]
            self.mean_molecular_wt += self.monomer_fraction[mo_type]*self.molecular_wt[mo_type]
        
        self.system_vol = self.total_sites*self.mean_molecular_vol
        #
        #self.population_density = (1E-24)*density/(mass_factor*self.mean_molecular_wt)
        
        self.total_effective_units = int(self.occupancy*self.total_sites)
        
        #self.total_effective_units = int(self.system_vol*self.population_density)
        
    def get_volumes(self):
        return self.mean_molecular_vol, self.system_vol, self.total_effective_units, self.mean_molecular_wt