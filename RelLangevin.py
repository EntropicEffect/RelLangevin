# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 15:06:12 2013

@author: Wiliam
"""
import tables
from pylab import *
import tests2

def En(p,m):
    return sqrt(m+ p.dot(p))

kb = 1 #GeV
c= 1   #Gev

dim =3
T = 13. #Gev
m = 1.3 #Gev
t = 50. #Gev-1
N = 0.5e4
dt = t/N #Gev-1


D = (T**2)/m

NumParticles = 100000

table = tables.openFile('NP100kT13t50','w')
atom = tables.Atom.from_dtype(dtype('Float64'))
pstore = table.createCArray(table.root,'Momentum',atom,[N,NumParticles,dim])

p = ones([2,NumParticles,dim])*2
pstore[0] = p[0]
#x= zeros([NumParticles,N,dim])

sqrtD = sqrt(D*dt)
m2 = m**2


for i in range(int(N-1)):
    dW = sqrtD*randn(1,NumParticles,dim)
    
    for n in range(NumParticles):
        p[1][n] = p[0][n] - D*p[0][n]*dt/(2*En(p[0][n],m2)*T) + dW[0][n]
      # p[1][n] = p[0][n] - Lorentz(p[0][n],m)*D*dt*p[0][n] + dW[0][n] 
        
        #x[n,i] = p[1][n]/(sqrt(p[1][n].dot(p[1][n]) + m2))
        
    pstore[i+1] = p[1]
    p[0] = p[1]
    
#tests2.MED2(p[0],T/10,10*T,T,m,NumParticles,dim)    
table.close()

