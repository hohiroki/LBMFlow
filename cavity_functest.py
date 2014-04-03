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
#print 'w'+str(w)
print 'u.shape:'+str(u.shape)
#print 'u'+str(u)
print 'c.shape:'+str(c.shape)
#print 'c'+str(c)
print 'rho.shape:'+str(rho.shape)
#print 'rho:'
#print rho
#print 'weights:'+str(w)
print 'u0:'+str(u0)

np.random.seed(1)
u = np.random.random((2,M,N))

print 'u changed to random:'
print u.shape
#print u

np.random.seed(1)
rho = np.random.random((M,N))

print 'rho changed to random:'
print rho.shape
#print rho

feq1 = otherEquilibrium(rho,u,w,c)
feq2 = otherEquilibrium2(rho,u,w,c)

#print 'feq1:'
#print feq1

#print 'feq2:'
#print feq2

feqdiff = feq1-feq2

print 'diff feq1-feq2'
#print feqdiff
print 'sum feqdiff'
print feqdiff.sum()

print '----'

# vel1 = velocity(f,rho,c)
# vel2 = otherVelocity(f,rho,c)
#
# print 'veldiff:'
# print (vel1-vel2).sum()

print

np.random.seed(1)
f = np.random.random((9,M,N))

print 'f changed to random'
print f.shape

np.random.seed(1)
feq = np.random.random((9,M,N))

print 'feq changed to random'
print feq.shape


# col1 = collision(f,feq,fnext,omega)
# col2 = collisionTwo(f,feq,fnext,omega)
# #col3 = otherCollision(f,feq,omega)
# print 'coldiff:'
# print (col1-col2).sum()

# den1 = density(f)
# den2 = densityTwo(f)
#
# print 'dendiff:'
# print (den1-den2).sum()

# stream1 = streaming(f)
# stream2 = streamingTwo(f)
# print 'streamdiff:'
# print (stream1-stream2).sum()


