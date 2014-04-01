from __future__ import division

__author__ = 'mh'

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
#import matplotlib as mp

from pylab import *

from cavity_params import *
from cavity import *

nsteps = 20000

rhosum_bound = np.zeros((5,nsteps+1))     # 0-4 Tot E N W S
rhosum_sect = np.zeros((2,nsteps+1))      # 0-1 horz vert

rhosum_bound[0,0] = rho.sum()           # total
rhosum_bound[1,0] = rho[:,-1].sum()     # E
rhosum_bound[2,0] = rho[0,:].sum()      # N
rhosum_bound[3,0] = rho[:,0].sum()      # W
rhosum_bound[4,0] = rho[-1,:].sum()     # S

rhosum_sect[0,0] = rho[M/2,:].sum()   # horz
rhosum_sect[1,0] = rho[:,N/2].sum()   # vert

steparray = np.linspace(0,nsteps+1,num=nsteps+1)

print 'steparray.shape'
print steparray.shape

print 'rhosum_bound.shape'
print rhosum_bound.shape
print rhosum_bound[0].shape

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
    u = fasterVelocity(f,rho,c)

    rhosum_bound[0,step+1] = rho.sum()           # total
    rhosum_bound[1,step+1] = rho[:,-1].sum()     # E
    rhosum_bound[2,step+1] = rho[0,:].sum()      # N
    rhosum_bound[3,step+1] = rho[:,0].sum()      # W
    rhosum_bound[4,step+1] = rho[-1,:].sum()     # S

    rhosum_sect[0,step+1] = rho[M/2,:].sum()   # horz
    rhosum_sect[1,step+1] = rho[:,N/2].sum()   # vert

    if (step % 1000) == 0:
        #print 'sum f:'+str(f.sum())
        print 'sum rho:'+str(rho.sum())
        #print 'sum feq:'+str(feq.sum())
    if (step % 1000) == 0:
        print 'Done step:'+str(step)


##########################################################

plt.figure(7)
Q = quiver( np.flipud(u[0]),np.flipud(u[1]))
# qk = quiverkey(Q, 0.5, 0.92, 2, r'$2 \frac{m}{s}$', labelpos='W',
#                fontproperties={'weight': 'bold'})
l,r,b,t = axis()
dx, dy = r-l, t-b
axis([l-0.05*dx, r+0.05*dx, b-0.05*dy, t+0.05*dy])

title('Velocity field')
plt.show()

# get the norm of the velocity
norm_u = np.sqrt(u[0]**2 + u[1]**2)

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.imshow(norm_u, cmap=cm.Reds)
plt.show()


plt.figure(3)

plt.plot(steparray,rhosum_bound[1],alpha = 0.5, linewidth=2., color='b')
plt.plot(steparray,rhosum_bound[2],alpha = 0.5, linewidth=2., color='r')
plt.plot(steparray,rhosum_bound[3],alpha = 0.5, linewidth=2., color='b')
plt.plot(steparray,rhosum_bound[4],alpha = 0.5, linewidth=2., color='b')

plt.plot(steparray,rhosum_sect[0], 'k--', alpha = 0.9, linewidth=3.)#color='k--')
plt.plot(steparray,rhosum_sect[1],'g:', alpha = 0.5, linewidth=3.)#color='g:')

#plt.plot(steparray,rhosum_bound[0], alpha = 0.5, linewidth=3., color='c')

plt.title('rho sum at boundaries E=b N=r W=b S=b and horz=k vert=g tot=c')
plt.show()


plt.figure(6)
plt.plot(steparray[:N],norm_u[M/2,:], alpha = 0.5, linewidth=3.)
plt.plot(steparray[:M],norm_u[:,N/2], alpha = 0.5, linewidth=3.)
plt.title('veocity profile horz=blue vert=green')
plt.show()