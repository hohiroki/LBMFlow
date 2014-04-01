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
print 'c'+str(c)
print 'rho.shape:'+str(rho.shape)
#print 'rho:'
#print rho
print 'weights:'+str(w)
print 'u0:'+str(u0)

nsteps = 2000


print 'sum f:'+str(f.sum())
print 'sum rho:'+str(rho.sum())
#print 'sum feq:'+str(feq.sum())

for step in range(0,nsteps):

    feq = equilibrium(rho,u,w,c)
    f = collision(f,feq,omega)

    f = streaming(f)

    f = westBoundary(f)
    f = eastBoundary(f)
    f = southBoundary(f)
    f,rhoN = northBoundary(f,u0)

    rho = density(f,rhoN)
    u = velocity(f,rho,c)

    if (step % 100) == 0:
        print 'sum f:'+str(f.sum())
        print 'sum rho:'+str(rho.sum())
        print 'sum feq:'+str(feq.sum())
    if (step % 100) == 0:
        print 'Done step:'+str(step)


# u = np.zeros((2,M,N))
# u[0,:,50]=np.arange(0.,1.,0.01)
# u[1,:,50]=np.arange(0.,1.,0.01)
# u[0,0,:]=np.arange(0.,1.,0.01)
# u[1,0,:]=np.arange(0.,1.,0.01)



fig = plt.figure(1)
ax = fig.add_subplot(111)
# ax.imshow(u[0,:,:]/np.linalg.norm(u), cmap='gray', interpolation='nearest')
# ax.imshow(u[1,:,:]/np.linalg.norm(u), cmap='gray', interpolation='nearest')
# ax.imshow(u[0,:,:], cmap='gray', interpolation='nearest')
# ax.imshow(u[1,:,:], cmap='gray', interpolation='nearest')
ax.imshow(np.sqrt(u[0]**2 + u[1]**2), cmap=cm.Reds)
plt.show()

plt.figure(2)
plt.streamplot(np.arange(0,M),np.arange(0,N),u[0,:,:],u[1,:,:])
# plt.plot(u[0,:,:])
# plt.plot(u[1,:,:])

# plt.figure()
# plt.contour(u[0,:,:])
# plt.contour(u[1,:,:])

plt.show()