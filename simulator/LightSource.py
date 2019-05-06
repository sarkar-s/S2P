# Object defining the light properties
import numpy as np
import math
import sys

class LightSource(object):
    def __init__(self,input_dict):
        self.lightmode = input_dict['Lighting mode (Linear ramp or Pulse)']
        self.KI = float(input_dict['Excitation rate constant'])
        self.time_limit = float(input_dict['Curing time (in seconds)'])

        self.reaction_type = input_dict['Initiation type (unimolecular or bimolecular)']

        if self.reaction_type == 'bimolecular':
            self.k_fl = 1.0/float(input_dict['Excitation lifetime of the initiator'])
            vol_1 = float(input_dict['Initiator volume (in A^3)'])
            vol_2 = float(input_dict['Co-initiator volume (in A^3)'])
            
            #area_1 = float(input_dict['Initiator surface area (in A^2)'])
            #area_2 = float(input_dict['Initiator surface area (in A^2)'])
            
            self.radius_1 = (vol_1*3.0/(math.pi*4.0))**(1.0/3.0)
            self.radius_2 = (vol_2*3.0/(math.pi*4.0))**(1.0/3.0)

            # Center to center distance
            self.sigma12 = self.radius_1 + self.radius_2

            self.molWt_1 = float(input_dict['Initiator molecular weight'])
            self.molWt_2 = float(input_dict['Co-initiator molecular weight'])

            self.reduced_mass = self.molWt_1*self.molWt_2/(self.molWt_1+self.molWt_2)
            self.reduced_mass *= (1E-26)/6.022

            self.temp = float(input_dict['Temperature (in K)'])

            self.eta = float(input_dict['Viscosity (in Pa.s)'])
            self.k_fl *= 0.1/self.eta
            # kB - Boltzmann factor
            self.kB = 1.38065E-23

            # combined velocity in Angstrom
            self.vel12 = (1E10)*np.sqrt((8.0*self.kB*self.temp)/(math.pi*self.reduced_mass))

            # Collision factor
            self.collision_factor = self.sigma12*self.vel12

            # Initial diffusivities 
            self.diffusivity_I = self.kB*self.temp/(6.0*math.pi)
            self.diffusivity_I *= (1.0/(self.eta*self.radius_1*(1E-10)))*(10**20)
            self.diffusivity_co_I = self.kB*self.temp/(6.0*math.pi)
            self.diffusivity_co_I *= (1.0/(self.eta*self.radius_2*(1E-10)))*(10**20)
            
        if self.lightmode == 'Linear ramp':
            self.ramp_time = float(input_dict['Ramp time (in seconds)'])

            self.initial_KI = 0.01*self.KI
            if self.ramp_time > 0:
                self.slope = 0.99*self.KI/self.ramp_time 
                self.nI_const = math.exp(-(self.initial_KI+0.5*self.slope*self.ramp_time)*self.ramp_time)      
        elif self.lightmode == 'Pulse':
            KI_string = input_dict['Pulse rate constants']
            KI_values = KI_string.split(',')
            
            time_string = input_dict['Pulse times (in seconds)']
            time_values = time_string.split(',')
            
            self.pulse_KIs = []
            self.pulse_times = []
                            
            self.pulse_set = range(0,len(KI_values))
            
            for d_idx in self.pulse_set:
                self.pulse_KIs.append(float(KI_values[d_idx]))
                self.pulse_times.append(float(time_values[d_idx]))

    def set_population(self,nI,n_co_I,system_vol):
        if self.reaction_type=='bimolecular':
            self.nI = nI
            self.n_co_I = n_co_I
            self.nI_excited = 0
            self.system_vol = system_vol
        elif self.reaction_type=='unimolecular':
            self.nI = nI
    
    def compute_total_diffusivity(self,d_factor_I,d_factor_co_I):
        self.D_12 = d_factor_I*self.diffusivity_I + d_factor_co_I*self.diffusivity_co_I
            
    def get_initiation_propensity(self):
        if self.reaction_type == 'unimolecular':
            propensity = self.lightfunction(time)*self.nI
        elif self.reaction_type == 'bimolecular':
            propensity = 4*math.pi*self.D_12*self.vel12*(self.sigma12**2)
            propensity *= 1.0/(self.system_vol*(4.0*self.D_12 + self.sigma12*self.vel12))
            propensity *= self.nI_excited*self.n_co_I
        return propensity
        
    def update_population(self):
        if self.reaction_type == 'unimolecular':
            self.nI += -1   
        elif self.reaction_type == 'bimolecular':
            self.nI_excited += -1
            self.n_co_I += -1

    def update_excitation(self):
        self.nI += -1
        self.nI_excited += 1
    
    def update_decay(self):
        self.nI += 1
        self.nI_excited += -1

    def get_excitation_propensity(self,time):
        propensity = self.lightfunction(time)*self.nI
        
        return propensity

    def get_decay_propensity(self):
        propensity = self.k_fl*self.nI_excited
                
        return propensity
        
    def lightfunction(self,time):
        
        if time > self.time_limit:
            this_KI = 0.0
        else:
            if self.lightmode == 'Linear ramp':
                # Linear ramp
                if time < self.ramp_time:
                    this_KI = self.initial_KI + self.slope*time
                else:
                    this_KI = self.KI
            elif self.lightmode == 'Instant':
                this_KI = self.KI
    
            elif self.lightmode == 'Pulse':
                time_low = 0.0
                
                if time>=self.time_limit:
                    this_KI = 0.0
                elif time==0.0:
                    this_KI = self.pulse_KIs[0]
                else:
                    for set_no in self.pulse_set[1:]:
                        if time>time_low and time<=self.pulse_times[set_no]:
                            this_KI = self.pulse_KIs[set_no-1]
                            break
                        
                        time_low = self.pulse_times[set_no]
        
        return this_KI

    def initiationFactor(self,current_time,last_time):
        if last_time <= self.ramp_time:
            nI_fac_low = math.exp(-(self.initial_KI+0.5*self.slope*last_time)*last_time)
        else:
            nI_fac_low = self.nI_const*math.exp(-self.KI0*(last_time-self.ramp_time))
            
        if current_time <= self.ramp_time:
            nI_fac_up = math.exp(-(self.initial_KI+0.5*self.slope*current_time)*current_time)
        else:
            nI_fac_up = self.nI_const*math.exp(-self.KI0*(current_time-self.ramp_time))
            
        return (nI_fac_low - nI_fac_up)