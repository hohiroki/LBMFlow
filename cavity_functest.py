from __future__ import division

__author__ = 'mh'

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from cavity_params import *
from cavity import *

print 'Reynolds:'+str(reynolds)
print 'Aspect:'+str(aspect)
print 'omega:'+str(omega)
print 'M:'+str(M)
print 'N:'+str(N)
print 'w.shape:'+str(w.shape)
print 'w'+str(w)
print 'u.shape:'+str(u.shape)
print 'u'+str(u)
print 'c.shape:'+str(c.shape)
print 'c'+str(c)
print 'rho.shape:'+str(rho.shape)
print 'rho:'
print rho
print 'weights:'+str(w)
print 'u0:'+str(u0)

u[:] = np.random.random((2,M,N))

print 'u changed to:'
print u.shape
print u

rho[:] = np.random.random((M,N))

print 'rho changed to:'
print rho.shape
print rho

feq1 = equilibrium(rho,u,w,c)
feq2 = otherEquilibrium(rho,u,w,c)

print 'feq1:'
#print feq1

print 'feq2:'
#print feq2

feqdiff = feq1-feq2

print 'diff feq1-feq2'
#print feqdiff
print 'sum diff'
print feqdiff.sum()