# Class for performing thermal motion of polymer chains
import math
import sys
import numpy as np

class ThermalMotion(object):
    def __init__(self,temp,sysVol,_bond_energy,_bond_angle_energy,_lattice_size,_lattice_h):
        # kB - Boltzmann factor
        self.kB = 1.38065E-23

        self.temp = temp
        self.sysVol = sysVol
        self.lattice_size = _lattice_size
        self.d_rms = self.lattice_size**2
        self.lattice_h = _lattice_h
        self.bond_energy = _bond_energy
        self.bond_angle_energy = _bond_angle_energy
        
        self.beta = self.kB*self.temp/(1.602E-19)
        
        # Energy coefficient in terms of kBT
        self.energy_eps = _bond_energy

        # Initial persistence length
        self.a0 = self.sysVol**(1.0/3.0)
        
        # Maximum radius size for FENE model in lattice unit
        self.r_max = (11.0)**(0.5)
        # Maximum bond stretch
        self.r_bond = (10.0)**(0.5)

        self.backUpTime = 0.0

        # Maximum FENE energy (for a junction point with functionality 2)
        #self.max_bond_energy = -0.5*67.5*self.energy_eps*math.log(1.0 - (self.r_bond/self.r_max)**2)
        self.max_bond_energy = self.energy_eps
        
        # Smallest value of energy change possible
        self.energy_quanta = self.energy_eps
        self.epr_quanta = self.energy_eps/self.beta
        
        #self.min_bond_energy = -0.5*67.5*self.energy_eps*math.log(1.0 - (2.0/self.r_max)**2)
        
        # Compute minimum energy bond angle
        self.minimum_energy_angle = math.acos(-1.0)
        
        self.energy_factor = math.exp(-self.energy_eps/self.beta)
        
        self.create_energy_change_dict()
        
        self.solvent_prob_keys = {}
    
    def compute_bond_angle_energies(self,_bond_angles):
        self.bond_angles = _bond_angles
        
        self.bond_angle_energies = {}

        for this_angle in self.bond_angles:
            if math.cos(this_angle)!=-1.0:
                self.bond_angle_energies[math.cos(this_angle)] = self.bond_angle_energy
            
        self.max_angle_energy = max(self.bond_angle_energies.values())
        
        self.create_maximum_energy_list()
        #self.create_MB_equilibrium_densities()
        
    def create_energy_change_dict(self):
        self.bond_energy_quanta = max(0.5*self.bond_energy,1e-8)
        self.bond_angle_energy_quanta = max(0.5*self.bond_angle_energy,1e-8)
        
        self.length_energy_factor = math.exp(-self.bond_energy/self.beta)
        self.angle_energy_factor = math.exp(-self.bond_angle_energy/self.beta)
        
        self.transition_odds = {}
        
        # Maximum energy change can happen when a crosslinked node is moved
        # this can change upto 4 bonds and 6 bond angles
        max_bond_idx = 8
        max_angle_idx = 12
        
        for i in xrange(-max_bond_idx,max_bond_idx+1):
            for j in xrange(-max_angle_idx,max_angle_idx+1):
                this_key = (i,j)
                energy_change = i*self.bond_energy_quanta + j*self.bond_angle_energy_quanta
                
                if energy_change>0:
                    self.transition_odds[this_key] = math.exp(-energy_change/self.beta)
                else:
                    self.transition_odds[this_key] = 1.0
        
    def create_maximum_energy_list(self):
        # This function stores the maximum possible energy values for:
        # (1) A free-end node in the network
        # (2) A node connected with 2 other nodes
        # (3) An active radical node connected at a junction with 3 other nodes
        # (3) A node connected to 4 other nodes or a junction point
        
        # The maximum amount of energy each node can contain
        self.maximum_energies = [0,0,0,0]
        
        # The maximum amount of energy change each node shift can cause
        self.max_displacement_energies = [0,0,0,0]
        self.max_index = [0,0,0,0]
        
        # Multipliers to get the transition probability
        self.transition_factors = [0,0,0,0]
        
        # Minimum escape probabilities
        self.minimum_escape_probs = [0,0,0,0]
        
        # Maximum energies of each type of network node
        # A free-end node can effect 1 bond length
        self.maximum_energies[0] = 0.5*self.max_bond_energy
        self.max_transition_factor = math.exp(-self.max_bond_energy/self.beta)
        # A middle-node can effect 2 bond lengths and 1 bond angles
        self.maximum_energies[1] = 0.5*2.0*self.max_bond_energy + 1.0*self.max_angle_energy
        # A Junction node with radical end can effect 3 bond lengths and 1 bond angles
        self.maximum_energies[2] = 0.5*3.0*self.max_bond_energy + 1.0*self.max_angle_energy
        # A Junction node can effect 4 bond lengths and 2 bond angles
        self.maximum_energies[3] = 0.5*4.0*self.max_bond_energy + 2.0*self.max_angle_energy

        # Maximum displacement energies of each type of network node
        # A free-end node can effect 1 bond length and 1 angle
        self.max_displacement_energies[0] = self.max_bond_energy + self.max_angle_energy
        # A middle-node can effect 2 bond lengths and 3 bond angles
        self.max_displacement_energies[1] = 2.0*self.max_bond_energy + 3.0*self.max_angle_energy
        # A Junction node with radical end can effect 3 bond lengths and 4 bond angles
        self.max_displacement_energies[2] = 3.0*self.max_bond_energy + 4.0*self.max_angle_energy
        # A Junction node can effect 4 bond lengths and 6 bond angles
        self.max_displacement_energies[3] = 4.0*self.max_bond_energy + 6.0*self.max_angle_energy
        
        # Transition factors to decide escape probability
        energy_factor = math.exp(-self.energy_quanta/self.beta)
        # A free-end node can effect 1 bond length and 1 angles
        if self.energy_quanta!=0.0:
            self.max_index[0] = int(self.max_displacement_energies[0]/self.energy_quanta)
            self.transition_factors[0] = (1.0-energy_factor**self.max_index[0])/(1.0-energy_factor)
        else:
            self.max_index[0] = 0.0
            self.transition_factors[0] = self.max_index[0]
            
        # A middle-node can effect 2 bond lengths and 3 bond angles
        if self.energy_quanta!=0.0:
            self.max_index[1] = int(self.max_displacement_energies[1]/self.energy_quanta)
            self.transition_factors[1] = (1.0-energy_factor**self.max_index[1])/(1.0-energy_factor)
        else:
            self.max_index[1] = 0.0
            self.transition_factors[0] = self.max_index[1]
            
        # A Junction node with radical end can effect 3 bond lengths and 4 bond angles
        if self.energy_quanta!=0.0:
            self.max_index[2] = int(self.max_displacement_energies[2]/self.energy_quanta)
            self.transition_factors[2] = (1.0-energy_factor**self.max_index[2])/(1.0-energy_factor)
        else:
            self.max_index[2] = 0.0
            self.transition_factors[2] = self.max_index[2]
            
        # A Junction node can effect 4 bond lengths and 6 bond angles
        if self.energy_quanta!=0.0:
            self.max_index[3] = int(self.max_displacement_energies[3]/self.energy_quanta)
            self.transition_factors[3] = (1.0-energy_factor**self.max_index[3])/(1.0-energy_factor)
        else:
            self.max_index[3] = 0.0
            self.transition_factors[3] = self.max_index[3]

        #self.compute_escape_probabilities()
        #self.compute_single_transition_odds()
        #self.create_transition_odds_dict()
                
    def compute_escape_probabilities(self):
        """ This function computes the probability of a node escaping
        from its current energy state to any other energy state
        """
        
        self.escape_probabilities = {}
        
        for node_type in range(0,4):
            self.escape_probabilities[node_type] = np.zeros(shape=(self.max_index[node_type]+1))
            idx_set = range(0,self.max_index[node_type]+1)
            
            for energy_idx in idx_set:
                this_energy = self.energy_quanta*energy_idx
                
                total_odds = 0.0
                #self.escape_probabilities[node_type][energy_idx] = 1.0 - math.exp(-energy_idx*self.energy_quanta/self.beta)/self.transition_factors[node_type]

                for inner_idx in idx_set:
                    final_energy = self.energy_quanta*inner_idx
                    total_odds += self.get_transition_prob(this_energy,final_energy)
                
                self.escape_probabilities[node_type][energy_idx] = total_odds/float(self.max_index[node_type]+1)
                
    def create_transition_odds_dict(self):
        self.transition_odds = {}
        self.transition_odds[0] = 1.0

        total_energy_idx = 2*int(self.max_displacement_energies[3]/self.energy_quanta)
        
        for idx in xrange(0,total_energy_idx+1):
            self.transition_odds[-idx] = 1.0
            self.transition_odds[idx] = math.exp(-0.5*idx*self.energy_quanta/self.beta)
        
    def get_transition_prob(self,current_E,final_E):
        energy_factor = math.exp(-(final_E - current_E)/self.beta)
        
        return min(1.0,energy_factor)

    def compute_single_transition_odds(self):
        """ There can be three types of events
        (1) Increase by an energy quanta
        (2) Decrease by an energy quanta
        (3) Energy remaining constant
        """
        sum_odds = math.exp(self.energy_quanta*self.beta)+1.0+math.exp(-self.energy_quanta*self.beta)
        
        # Probability of a reduction step
        self.reduction_odds = math.exp(self.energy_quanta*self.beta)/sum_odds
        
        # Probability of an increase step
        self.increase_odds = math.exp(-self.energy_quanta*self.beta)/sum_odds
        
        #print self.reduction_odds,self.increase_odds
        #sys.exit()
        
    def compute_equilibriation_steps(self,total_energy):
        N = int(total_energy/(self.energy_quanta*(self.reduction_odds - self.increase_odds)))
        
        return N

    def maximum_fluctuation_propensity(self,diffusivity):
        max_value = 6.0*diffusivity*max(self.lattice_propensity)
        
        return max_value

    def create_MB_equilibrium_densities(self):
        jn_MB_sum = (1.0 - math.exp(-13.0*0.5*self.energy_eps/self.beta))/(1.0 - math.exp(-0.5*self.energy_eps/self.beta))
        self.jn_MB_probs = np.zeros(shape=(13))
        for idx in range(0,13):
            self.jn_MB_probs[idx] = math.exp(-0.5*float(idx)*self.energy_eps/self.beta)/jn_MB_sum
            
        nn_MB_sum = (1.0 - math.exp(-5.0*0.5*self.energy_eps/self.beta))/(1.0 - math.exp(-0.5*self.energy_eps/self.beta))
        self.nn_MB_probs = np.zeros(shape=(13))
        for idx in range(0,5):
            self.nn_MB_probs[idx] = math.exp(-0.5*float(idx)*self.energy_eps/self.beta)/nn_MB_sum
            
        nr_MB_sum = (1.0 - math.exp(-2.0*0.5*self.energy_eps/self.beta))/(1.0 - math.exp(-0.5*self.energy_eps/self.beta))
        self.nr_MB_probs = np.zeros(shape=(13))
        for idx in range(0,2):
            self.nr_MB_probs[idx] = math.exp(-0.5*float(idx)*self.energy_eps/self.beta)/nr_MB_sum
            
        #print self.nr_MB_probs, self.nn_MB_probs, self.jn_MB_probs
        #sys.exit()
        
    def compute_monomer_self_diffusion(self,unit_vol,mol_vols,eta):
        self.lattice_h = 0.5*(unit_vol**(1.0/3.0))

        # Diffusion coefficient in A^2 s^-1
        self.Diff = []
        self.lattice_propensity = []
        
        for this_vol in mol_vols:
            # Spherical volume size
            r_s = ((3.0/(4.0*math.pi))*this_vol)**(1.0/3.0)
            r_s *= 1E-10
            
            if eta == "infinite":
                self.Diff.append(0.0)
                self.lattice_propensity.append(0.0)
            else:
                eta = float(eta)
                self.Diff.append(self.kB*self.temp/(6.0*math.pi))
                self.Diff[-1] *= (1E+20/(float(eta)*r_s))

                self.lattice_propensity.append(self.Diff[-1]/(self.lattice_h**2))
            
            self.mean_lattice_propensity = sum(self.lattice_propensity)/float(len(mol_vols))
            
    def get_monomer_self_diffusion(self,hopping_rate):
        # Diffusion coefficient in A^2 s^-1
        self.lattice_propensity = []
        
        for hop in hopping_rate:
            self.lattice_propensity.append(hop)
            
        self.mean_lattice_propensity = sum(self.lattice_propensity)/float(len(hopping_rate))
                
    def compute_mean_diffusivity(self,mol_fraction):
        self.mean_diff = 0.0
        self.mean_propensity = 0.0
        
        for mo_type in xrange(0,len(mol_fraction)):
            self.mean_diff += mol_fraction[mo_type]*self.Diff[mo_type]
            self.mean_propensity += mol_fraction[mo_type]*self.lattice_propensity[mo_type]
        
    def get_lattice_propensity(self,surface_site_count):
        if surface_site_count == 0:
            surface_site_count = 26.0
            
        molecular_weight_factor = surface_site_count/26.0

        new_propensity = self.lattice_propensity/molecular_weight_factor

        return new_propensity
        
    def get_no_of_points(self,timeStep):
        self.get_melt_reptation_time()

        timeStep += self.backUpTime

        no_of_moves = int(float(self.total_chains)*(timeStep/self.solTime))
        
        no_of_moves = int(timeStep/self.solTime)
        
        if no_of_moves > 0:
            self.backUpTime = 0.0
        else:
            self.backUpTime = timeStep
            
        print 'Moves : ', no_of_moves

        return no_of_moves
        
    def get_chain_density(self,chains,cross_links,used_monomers):
        if cross_links == 0:
            self.total_chains = chains
        else:
            self.total_chains = self.funcNum*cross_links

        self.chain_density = self.total_chains/self.sysVol
        
        # Average degree of polymerization
        self.mean_polymerization = used_monomers/self.total_chains

        # Compute persistence length
        self.get_persistence_length()

        # Mean arc length
        self.mean_l = self.mean_polymerization*(1.0/self.a)

        # Compute disengagement time
        self.T_d = (self.mean_l**2)/(self.Diff*math.pi**2)
    
    def get_persistence_length(self):
        self.a = self.a0*(self.total_chains)**(-3.0/5.0)

        print 'Persistence length : ', self.a

    # Compute the time to next reptation event 
    def get_melt_reptation_time(self):
        self.d_rms = self.lattice_h**2
        
        # First approximation
        t_1 = self.d_rms*(self.mean_l/(2*self.a*self.Diff))

        # Find the time for reptation
        self.solTime = brentq(self.displacement_function,t_1,self.T_d)
        
        self.solTime *= (1.0/self.total_chains)
        
        return self.solTime


    # The propensity to have a thermal event within time interval dt
    def get_thermal_propensity(self,movable_nodes):
        # Propensity for one polymer unit to move by one lattice distance
        # in a given time interval
        if movable_nodes == 0:
            movable_nodes = 1

        propensity = 6.0*self.Diff/(self.d_rms*movable_nodes)

        # This is multiplied by the no. of units in the network
        propensity *= movable_nodes

        return propensity
    
    def get_FENE_energy_change(self,oldbox,newbox,refbox):
        # Old distance
        r_old, r_new = 0.0, 0.0
        
        for i in range(0,3):
            short1 = abs(refbox[i] - oldbox[i])
            short2 = abs(abs(refbox[i] - oldbox[i]) - self.lattice_size)
            if short1>short2:
                short1 = short2
                
            r_old += short1**2

            short1 = abs(refbox[i] - newbox[i])
            short2 = abs(abs(refbox[i] - newbox[i]) - self.lattice_size)
            if short1>short2:
                short1 = short2

            r_new += short1**2

        r_old = r_old**(0.5)
        # Initial energy        
        U_initial = -0.5*67.5*self.energy_eps*math.log(1.0 - (r_old/self.r_max)**2)

        r_new = r_new**(0.5)    
        U_final = -0.5*67.5*self.energy_eps*math.log(1.0 - (r_new/self.r_max)**2)
            
        U_diff = U_final - U_initial
        
        return U_diff 

    def get_FENE_energy(self,idx,moveIdx):
        r_value = 0.0
        
        for i in range(0,3):
            short1 = abs(idx[i] - moveIdx[i])
            short2 = abs(abs(idx[i] - moveIdx[i]) - self.lattice_size)
            if short1>short2:
                short1 = short2
                
            r_value += short1**2

        r_value = (r_value)**(0.5)

        U_value = -0.5*67.5*self.energy_eps*math.log(1.0 - (r_value/self.r_max)**2)

        return U_value
    
    def get_bond_energy_change(self,oldbox,newbox,refbox):
        # Old distance
        r_old, r_new = 0.0, 0.0
        
        for i in range(0,3):
            short1 = abs(refbox[i] - oldbox[i])
            short2 = abs(abs(refbox[i] - oldbox[i]) - self.lattice_size)
            if short1>short2:
                short1 = short2
                
            r_old += short1**2

            short1 = abs(refbox[i] - newbox[i])
            short2 = abs(abs(refbox[i] - newbox[i]) - self.lattice_size)
            if short1>short2:
                short1 = short2

            r_new += short1**2

        # Initial energy
        if math.sqrt(r_old) != 3.0:
            U_initial = self.energy_eps
        else:
            U_initial = 0.0
        
        if math.sqrt(r_new) != 3.0:
            U_final = self.energy_eps
        else:
            U_final = 0.0
            
        U_diff = U_final - U_initial
        
        return U_diff, r_new, r_old
        
    def get_bond_energy(self,idx,nextIdx):
        r_value = 0.0
        
        for i in range(0,3):
            short1 = abs(idx[i] - nextIdx[i])
            short2 = abs(abs(idx[i] - nextIdx[i]) - self.lattice_size)
            if short1>short2:
                short1 = short2
                
            r_value += short1**2

        r_value = (r_value)**(0.5)
        
        if r_value != 3.0 and r_value != 0.0:
            U_value = self.energy_eps
        else:
            U_value = 0.0

        return U_value
    
    def get_angle(self,start_node,center_node,end_node):
        vector_1, vector_2 = [], []
        for i in range(0,3):
            vector_1.append(start_node[i]-center_node[i])
            vector_2.append(end_node[i]-center_node[i])
            
        d_prod = np.inner(vector_1,vector_2)
        lv1 = np.linalg.norm(vector_1)
        lv2 = np.linalg.norm(vector_2)
            
        cos_theta = d_prod/(lv1*lv2)

        if abs(cos_theta)>1.0:
            cos_theta += -np.sign(cos_theta)*np.finfo(float).eps
            
        return cos_theta
    
    def get_angle_energy(self,cos_bond_angle):
        angle_difference = cos_bond_angle + 1.0
        
        #energy = 0.5*self.bond_angle_energy*(angle_difference**2)
        if angle_difference!=0.0:
            energy = self.bond_angle_energy
        else:
            energy = 0.0
        
        return energy
    
    def get_angle_energy_change(self,start_node,end_node,old_box,new_box):
        angle_initial = self.get_angle(start_node,old_box,end_node)
        angle_final = self.get_angle(start_node,new_box,end_node)

        U_angle_diff = self.get_angle_energy(angle_final)-self.get_angle_energy(angle_initial)
        
        return U_angle_diff
    
    def energy_propensity_factor(self,current_energy,tr_idx):
        energy_factor = math.exp(-current_energy/self.beta)*self.transition_factors[tr_idx]
        #max_idx = int(self.max_displacement_energies[tr_idx]/(0.5*self.energy_eps))
        #current_idx = int(current_energy/(0.5*self.energy_eps))
        #
        #updated_idx = max_idx - current_idx
        #
        #if updated_idx<0 or updated_idx>max_idx:
        #    print 'Energy index error: ', updated_idx, current_idx, max_idx
        #    
        #energy_factor = math.exp(0.5*float(updated_idx)*self.energy_eps/self.beta)
        
        return energy_factor