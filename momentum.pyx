# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 23:38:24 2013

@author: Wiliam
"""

import numpy as np
cimport numpy as np
ctypedef np.float64_t dtype_t
cimport cython
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)

def pupdate(np.ndarray[dtype_t,ndim=3] p,int N,double mu,double T,double dt,int NumParticles,int dim,double m2):
    cdef double ener
    cdef int i,n,c,d
    d = 0
    c = 1
    for i in xrange(N-1):
        for n in xrange(NumParticles):
             ener  = np.sqrt(m2+ p[0][n].dot(p[0][n]))
             drag = - p[0][n]*dt*mu
             noise = np.random.randn(dim)
             #dp = drag + np.sqrt(2*mu*T*ener*dt/(1-T/ener))*noise
            # ener = np.sqrt(m2 + (p[0][n]+dp[:]).dot((p[0][n]+dp[:])))
             p[1][n] = p[0][n] + drag  + np.sqrt(2*mu*(ener+T)*T)*noise
        p[0] = p[1]
        d += 1
        if (d==N/10):
            d = 0
            print(c*10)
            c+=1
    return p[1]