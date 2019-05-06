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

class phonon_properties(object):
    def __init__(self,_omega, _g_omega, _MFP, _lambda_omega, _v_omega,_scale):
        # No. of bins in the histogram
        self.omega = _omega
        self.g_omega = _g_omega
        self.MFP = _MFP
        self.lambda_omega = _lambda_omega
        self.v_omega = _v_omega
        self.scale = _scale
        
        self.h_bar = 1.054E-34
        
        # kBT in eV
        self.kBT = 0.0259
        
        # cm to ev factor
        self.cm2ev = 1.24E-4
        
        # Cv and Kappa
        self.Cv, self.kappa, self.alpha = {}, {}, {}
        
        # Temperature range
        self.temp = np.linspace(1,500,500)
        
        for i in self.temp:
            self.Cv[i] = 0.0
            self.kappa[i] = 0.0
            self.alpha[i] = 0.0
        
        self.mean_v = 0.0
        self.mean_l = 0.0
        self.mean_v_l = 0.0
        
    def compute_cv_kappa(self):
        self.ratio = []
        
        total = 0.0
        
        for i in xrange(0,len(self.omega)-1):
            d_omega = self.omega[i+1]-self.omega[i]
                
            g_mean = 0.5*(self.g_omega[i]+self.g_omega[i+1])
                
            omega_mean = 0.5*(self.omega[i]+self.omega[i+1])
                
            true_g = 0.5*(self.g_omega[i]*self.omega[i]**2+self.g_omega[i+1]*self.omega[i+1]**2)
            
            self.mean_v += true_g*d_omega*self.v_omega[i]
            
            self.mean_l += true_g*d_omega*self.MFP[i]
            
            self.mean_v_l += true_g*d_omega*self.v_omega[i]*self.MFP[i]
            
            total += d_omega*true_g
            
        self.mean_v *= 1.0/total
        self.mean_l *= 1.0/total
        self.mean_v_l *= 1.0/total
        
        for t in self.temp:
            t_factor = self.kBT*t/300.0
        
            for i in xrange(0,len(self.omega)-1):
                d_omega = self.omega[i+1]-self.omega[i]
                
                g_mean = 0.5*(self.g_omega[i]+self.g_omega[i+1])
                
                omega_mean = 0.5*(self.omega[i]+self.omega[i+1])
                
                x = self.cm2ev*omega_mean/t_factor
                
                true_g = 0.5*(self.g_omega[i]*self.omega[i]**2+self.g_omega[i+1]*self.omega[i+1]**2)
                
                #total += d_omega*true_g
                
                self.ratio.append(x)
            
                factor = x*self.scale
                
                try:
                    integral = (factor**2)*math.exp(factor)/((math.exp(factor)-1)**2)
                except OverflowError:
                    print factor
                    sys.exit()
                    
                self.Cv[t] += d_omega*true_g*integral
                
                self.kappa[t] += (1.0/3.0)*d_omega*true_g*integral*self.MFP[i]*self.v_omega[i]

                self.alpha[t] += d_omega*true_g*self.MFP[i]*self.v_omega[i]/(math.exp(factor)-1)
                
            #self.kappa[t] = (1.0/3.0)*self.Cv[t]*self.mean_v_l#*self.mean_l

    def write_results(self):
        tD_idx = 0
        
        for i in xrange(0,len(self.temp)-1):
            tD_idx = i
            
            key1, key2 = self.temp[i], self.temp[i+1]
            
            if self.Cv[key1]<=0.9517 and self.Cv[key2]>0.9517:
                break
                
        factor = self.temp[tD_idx]
        
        print 'Debye Temperature: ', self.temp[tD_idx]
        
        ofile = open('phonon_properties.csv','w')
        
        print >> ofile, 'Temperature',',','Cv',',','kappa'
        
        for key in self.temp:
            print >> ofile, key/factor,',',str(self.Cv[key]),',',str(self.kappa[key]),',',str(self.alpha[key])
            
        ofile.close()
        
        ofile2 = open('modulus.csv','w')
        
        print >> ofile2, 'Modulus',',',self.mean_v**2
        
        ofile2.close()