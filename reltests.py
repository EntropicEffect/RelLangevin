# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 22:41:33 2013

@author: Wiliam
"""
from pylab import *
import scipy.special.kv as bessel
k =1.38e-23
def MED(p,m,c,T,numparticls):
    theta = (k*T/(m*c**2))
    c = 1/(theta*bessel(theta/2))
    p = p/m 
    gamma = zeros([1,numparticles])
    for i in range(numparticles):
        gamma[i] = 1/sqrt((1-(p[i].dot(p[i])/c**2)))
    gf =  100
    dg = 0.1
    GammaDist = arange(0,100,df)    
    Juttner = theta*GammaDist**2*(sqrt(1-(1-(1/GammaDist))))*e**(-GammaDist/theta)
    plot(GammaDist,Juttner)
        
        
    
    
    