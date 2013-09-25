# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 22:41:33 2013

@author: Wiliam
"""
from pylab import *
import scipy.special as sp
from scipy.special import gamma as Gamma
from scipy.special import hyp0f1 as hyperGeom
k =1

    
def MED2(p,T,m,n,dim):  
    """ p,T,m,n,dim """
    de = T/15
    ef = 11*T
    normconst = (2./m*k*T)**((dim-1)/2.)/sqrt(pi) *sp.kv((dim+1)/2.,(1.*m)/(k*T))*Gamma((dim/2.))
    energy = zeros([1,n])
    enrange = arange(0,ef,de)
    Juttner = zeros([1,len(enrange)])
    m2 = m**2
    dim = double(dim)
    for i in range(n):
        energy[0,i] = sqrt((p[i].dot(p[i]))+m2) 
    g= 0 
   
    for E in enrange:
        Juttner[0,g] = (n*de)*(1/normconst)*((E**2 -m**2)/(m**2))**(dim/2.)*(E/(E**2-m**2))*e**(-E/(k*T))
        g+= 1
    weights = ones([1,n])
    y,binEdges=np.histogram(energy[0],bins=len(enrange),range=(0,ef),weights=weights[0])
    menStd = sqrt((y-(y/n)))
    errorbar(enrange,y, yerr=menStd, fmt='ro') 
    plot(enrange,Juttner[0]) 
    numwithin = 0
    for i in range(len(enrange)):
        if (y[i]+menStd[i] >= Juttner[0,i]) and (y[i]-menStd[i] <= Juttner[0,i]):
            numwithin+= 1
    print((numwithin*100)/len(enrange), "% are within 1 standard dev")        
    
    
def boostDist(p,T,m,n,dim,beta):
    dp = 0.5
    pf = 60
    m2 = m**2
    gamma = 1/(sqrt(1-beta**2))
    energy = zeros([1,n])
    for i in range(n):
        energy[0,i] = sqrt(p[i].dot(p[i])+m2)
    for i in range(n):
        p[i,2] = gamma*p[i,2] - beta*gamma*sqrt(p[i].dot(p[i]) + m2)
    normconst = sqrt(pi)*(2./m*k*T)**((dim-1)/2.)/sqrt(pi) *sp.kv((dim+1)/2.,(1.*m)/(k*T))*Gamma((dim/2.))
    momentum = zeros([1,n])
    momrange = arange(-2*pf,pf,dp)
    gt = gamma/T
    Dist = zeros([1,n])
    dim = double(dim)
    for i in range(n):
        momentum[0,i] = p[i,2]
    g= 0
    for a in momrange:
        Dist[0,g] = (n*dp)*(1/normconst)*(hyperGeom((dim+2)/2,0.25*(gamma*beta*a/T)**2))*e**(-gt*sqrt(a**2 + m2))
        g += 1
    weights = ones([1,n])
    y,binEdges=np.histogram(momentum[0],bins=len(momrange),range=(-2*pf,pf),weights=weights[0])
    menStd = sqrt((y-(y/n)))
    errorbar(momrange,y, yerr=menStd, fmt='ro') 
    plot(momrange,Dist[0]) 
    
#*(1/hyperGeom((dim+2)/2,0.25*(gamma*beta*a/T)**2))
        
        
    
    
    