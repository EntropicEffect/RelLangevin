# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 22:41:33 2013

@author: Wiliam
"""
from pylab import *
import scipy.special as sp
from scipy.special import gamma as Gamma
from scipy.special import hyp0f1 as hyperGeom
from scipy.integrate import romberg as numIntegral
k =1
def JuttnerDist(x,m,dim,gamma,T,beta):
    m2 =m**2
    return (x**2 + m2)**(dim/4)*sp.kv(dim/2,sqrt(x**2 +m2)*gamma/T)*e**(gamma*beta[0]*x/T)
def BoostVector(beta,m,p):
    m2 = m**2
    N = int(len(beta))
    B = sqrt(beta.dot(beta))
    if B >= 1:
        print('Slow down, Einstein!')
    E = sqrt(p.dot(p) + m2)
 
    fourVector = zeros([1,N+1])
    fourVector[0,0] = E
    fourVector[0,1:] = p
    gamma = 1/sqrt(1-B**2)
    bmatrix = zeros([N+1,N+1])
    bmatrix[0,0] = gamma
    bmatrix[1:N+1,0] = -gamma*beta
    bmatrix[0,1:N+1] = -gamma*beta
    for i in range(1,N+1):
        bmatrix[i,i] =  1 + (gamma-1)*(beta[i-1]**2/B**2)
    
    for c in range(1,N+1):
        for r in range(1,N+1):
            if r == (c):
                continue
            bmatrix[r,c] = (gamma-1)*beta[c-1]*beta[r-1]/B**2


    return bmatrix
        
    
def MED2(p,T,m,n,dim):  
    """ p,T,m,n,dim """
    de = T/50 
    ef = 4.5
    normconst = (2./m*k*T)**((dim-1)/2.)/sqrt(pi) *sp.kv((dim+1)/2.,(1.*m)/(k*T))*Gamma((dim/2.))
    energy = zeros([1,n])
    enrange = arange(m+0.1,ef,de)
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
    y,binEdges=np.histogram(energy[0],bins=len(enrange),range=(m+0.1,ef),weights=weights[0])
    menStd = sqrt((y-(y/n)))
    errorbar(enrange,y, yerr=menStd, fmt='ro') 
    plot(enrange,Juttner[0]) 
    numwithin = 0
    for i in range(len(enrange)):
        if (y[i]+menStd[i] >= Juttner[0,i]) and (y[i]-menStd[i] <= Juttner[0,i]):
            numwithin+= 1
    print((numwithin*100)/len(enrange), "% are within 1 standard dev")        
    
    
def boostDist(p,T,m,n,dim,beta):
    dp = 0.2
    pf = 15
    m2 = m**2
    v = zeros(n)
    
    gamma = 1/(sqrt(1-beta[0]**2))
    for i in range(n):
#        fourVector = zeros([1,dim+1])
#        fourVector[0,0] = sqrt(p[i].dot(p[i]) + m2)
#        fourVector[0,1:dim+1] = p[i]
#
#        fourVector = BoostVector(beta,m,p[i]).dot(fourVector.T)
#    
#        for k in range(dim):
         v[i] = p[i][2]/sqrt(m2+ p[i].dot(p[i]))
         v[i] = (v[i] + beta[0])/(1+beta[0]*v[i])
         p[i][2] = m*v[i]/(sqrt(1-v[i]**2))
   
    normconst = 1/((2./m*k*T)**((dim-1)/2.)/sqrt(pi) *sp.kv((dim+1)/2.,(1.*m)/(k*T)))#*Gamma((dim/2.)))
    momentumx = zeros([1,n])
    momentumy = zeros([1,n])
    momentumz = zeros([1,n])
    momentum = zeros([1,n])
    momrange = arange(-pf,pf,dp)
    momlength = arange(0,pf,dp)
    gt = gamma/T
    Dist = zeros([1,len(momlength)])
    dim = double(dim)
    for i in range(n):
        momentumx[0,i] = p[i,0]
        momentumy[0,i] = p[i,1]
        momentumz[0,i] = p[i,2]
        momentum[0,i] =  sqrt(sum(p[i]**2))
    g= 0
    for a in momlength:
        Dist[0,g] = normconst*(n*dp)*a**2*hyperGeom((5/2),0.25*(gamma*beta[0]*a/T)**2)*e**(-gt*sqrt(a**2 + m2))
        g += 1
    weights = ones([1,n])
    subplot(311)
    y,binEdges=np.histogram(momentumx[0],bins=len(momrange),range=(-pf,pf),weights=weights[0])
    menStd = sqrt((y-(y/n)))
    errorbar(momrange,y, yerr=menStd, fmt='ro') 
    testplot = 800*e**(-gt*(sqrt(momrange**2 + m2) + momrange*beta[0]))
    
    plot(momrange,testplot)
    subplot(312)
    y,binEdges=np.histogram(momentumy[0],bins=len(momrange),range=(-pf,pf),weights=weights[0])
    menStd = sqrt((y-(y/n)))
    errorbar(momrange,y, yerr=menStd, fmt='ro')
    subplot(313)
    y,binEdges=np.histogram(momentumz[0],bins=len(momrange),range=(-pf,pf),weights=weights[0])
    menStd = sqrt((y-(y/n)))
    errorbar(momrange,y, yerr=menStd, fmt='ro')
   # plot(momrange,Dist[0]) 
    figure()
    y,binEdges=np.histogram(momentum[0],bins=len(momlength),range=(0,pf),weights=weights[0])
    menStd = sqrt((y-(y/n)))
    errorbar(momlength,y, yerr=menStd, fmt='ro')
  
    plot(momlength,Dist[0])
    show()
    
def pz(p,T,m,n,dim,beta):
    pcpy = copy(p)
    print(beta[0])
    pb = -10
    pf = 10
    dp = 0.1
    
    m2 = m**2
    gamma = 1/sqrt(1-beta[0]**2)
    norm = 1/numIntegral(JuttnerDist,-10,10,args=(m,dim,gamma,T,beta),divmax = 10)
    normconst = 1/(sqrt(2*pi*T)*m**((dim+1)/2)*gamma**((dim/2)-2)*sp.kv((dim+1)/2,m/T))
    momrange= arange(pb,pf,dp)
    print(norm)
    print(normconst)
    Juttner = zeros(len(momrange))
    Juttner2 = zeros(len(momrange))
    momz = zeros(n)
    v = zeros(n)
    for i in xrange(n):
         v[i] = pcpy[i][2]/sqrt(m2+ pcpy[i].dot(pcpy[i]))
         v[i] = (v[i] + beta[0])/(1+beta[0]*v[i])
         pcpy[i][2] = m*v[i]/(sqrt(1-v[i]**2))
   
    for i in xrange(n):
        momz[i] = pcpy[i][2]
    g = 0
    for x in momrange:
        Juttner[g] = n*dp*norm*(x**2 + m2)**(dim/4)*sp.kv(dim/2,sqrt(x**2 +m2)*gamma/T)*e**(gamma*beta[0]*x/T)
        Juttner2[g] = n*dp*normconst*(x**2 + m2)**(dim/4)*sp.kv(dim/2,sqrt(x**2 +m2)*gamma/T)*e**(gamma*beta[0]*x/T)
        g += 1
    weights = ones([1,n])
    y,binEdges=np.histogram(momz,bins=len(momrange),range=(pb,pf),weights=weights[0])
    menStd = sqrt((y-(y/n)))
    errorbar(momrange,y, yerr=menStd, fmt='ro') 
    plot(momrange,Juttner,hold=True)
    plot(momrange,Juttner2)
#*(1/hyperGeom((dim+2)/2,0.25*(gamma*beta*a/T)**2))
        
        
    
    
    
