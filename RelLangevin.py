# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 15:06:12 2013

@author: Wiliam
"""
import tables
from pylab import *
import tests2
from momentum import pupdate

kb = 1 #GeV
c= 1   #Gev
lamda = 12
sqrtL= sqrt(lamda)

dim =3
T = 0.45 #Gev
m = 1.5 #Gev
t = 15. #Gev-1
N = 1000
dt = t/N #Gev-
mu = (pi*sqrtL*T**2)/(2*m)
NumParticles = 2000

table = tables.openFile('NP10k','w')
atom = tables.Atom.from_dtype(dtype('Float64'))
pstore = table.createCArray(table.root,'Momentum',atom,[NumParticles,dim])

p = ones([2,NumParticles,dim])*randn(2,NumParticles,dim)*20

m2 = m**2


p = pupdate(p,N,mu,T,dt,NumParticles,dim,m2)   
#pupdate(p,int N,double sqrtD,double dt,int NumParticles,int dim,double m2)
#pstore[:] = p[1]
#tests2.MED2(p[0],T/10,10*T,T,m,NumParticles,dim)    
table.close()

