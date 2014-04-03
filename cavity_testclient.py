from __future__ import division

__author__ = 'mh'

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
#import matplotlib as mp
import scipy.io as sio

from pylab import *

from cavity_params import *
from cavity import *
import time
import datetime


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

    # gather statistics
    rhosum_bound[0,step+1] = rho.sum()           # total
    rhosum_bound[1,step+1] = rho[:,-1].sum()     # E
    rhosum_bound[2,step+1] = rho[0,:].sum()      # N
    rhosum_bound[3,step+1] = rho[:,0].sum()      # W
    rhosum_bound[4,step+1] = rho[-1,:].sum()     # S

    rhosum_sect[0,step+1] = rho[M/2,:].sum()   # horz
    rhosum_sect[1,step+1] = rho[:,N/2].sum()   # vert

        #print 'sum f:'+str(f.sum())
        #print 'sum rho:'+str(rho.sum())
        #print 'sum feq:'+str(feq.sum())
    if (step % 1000) == 0:
        print 'Done step:'+str(step)
        print 'sum rho:'+str(rho.sum())


##########################################################

# mat_dict = {'nsteps':np.asarray(nsteps),'u':u[0],'v':u[1],'u0':u0[0],'v0':u[1],
#             'rho0':np.asarray(rho0), 'reynolds':np.asarray(reynolds),
#             'visc':np.asarray(viscosity), 'vel':np.asarray(velocity),
#             'visc_lb':np.asarray(viscosity_lattice),
#             'vel_lb':np.asarray(velocity_lattice), 'rhosum_bound':rhosum_bound,
#             'rhosum_sect':rhosum_sect}

mat_dict = {'u':u[0],'v':u[1],
            'u0':u0[0],'v0':u[1],
            'nsteps':np.asarray(nsteps),
            'rho0':np.asarray(rho0),
            'reynolds':np.asarray(reynolds),
            'visc_ph':np.asarray(viscosity),
            'vel_ph':np.asarray(vel_ph),
            'visc_lb':np.asarray(viscosity_lattice),
            'vel_lb':np.asarray(velocity_lattice),
            'rhosum_bound':rhosum_bound,'rhosum_sect':rhosum_sect}

timestamp = time.time()
timestamp_formatted = datetime.datetime.fromtimestamp(timestamp).strftime('%Y-%m-%d_%H.%M.%S')

mat_filename = 'Clb_'+str(nsteps)+'_'+timestamp_formatted

sio.savemat(mat_filename,mat_dict,appendmat=True)

##########################################################

steparray = np.linspace(0,nsteps+1,num=nsteps+1)

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
title('norm(velocity)')
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