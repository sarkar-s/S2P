# Class to determine modulus of the polymer matrix

import os, sys
import numpy as np
import math
import cmath

class modulus(object):
    def __init__(self,temp,sysVol,func):
        # kB - Boltzmann factor
        self.kB = 1.38065E-23
        # Planck's constant
        self.h_bar = (6.626E-34)/(2.0*math.pi)
        
        # Temperature
        self.temp = float(input_dict['Temperature (in K)'])
        self.beta = 1.0/(self.kB*temp)

        # System volume
        self.sysVol = sysVol
        # Functionality of each cross-link
        self.funcPhi = func
        
        # High frequency moduli
        self.G_infi = float(input_dict['High frequency shear moduli (in GPa)'])
        self.K_infi = float(input_dict['High frequency bulk moduli (in GPa)'])
        self.E_infi = 9.0*self.K_infi*self.G_infi/(3.0*self.K_infi + self.G_infi)
        
        # Viscosity
        self.eta_0 = float(input_dict['Viscosity (in Pa.s)'])

        # Factor to modulus
        self.modFactor = self.kB*self.temp/self.sysVol        
                        
    def getGcGe(self,totalCL,CLperChain):
        self.Gc = (1E-6)*self.modFactor*totalCL*(1.0 - 2.0/self.funcPhi)

        self.Ge = (1E-6)*(4.0/7.0)*(self.modFactor*totalCL*CLperChain)

        return self.Gc, self.Ge
    
    def get_relaxation_time_constant(self,eta_factor):
        self.relax_tc = eta_factor*self.E_infi/self.eta_0