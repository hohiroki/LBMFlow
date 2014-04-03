from __future__ import division

__author__ = 'mh'

# D2Q9 Lattice Boltzmann

import numpy as np

# input data

# fluid physics
viscosity = 1.2*1e-3    # kinematic viscosity m2/s (oil at 15 deg C)

vel_ph = 6.0           # moving lid m/s

height = 0.2     # height of cavity m
width = 0.2      # width of cavity m

reynolds = vel_ph * height / viscosity    # macroscale reynolds number

# LBM values
aspect = width / height         # aspect ratio of domain

reynolds_lattice = reynolds     # keep same reynolds
viscosity_lattice = 0.01        # choose so we keep same reynolds
velocity_lattice = 0.1          # choose so we get 100x100 (here we ignore the actual velocity)
rho0 = 5.0

omega = 1./(3. * viscosity_lattice + 0.5)   # omega = 1/tau

N = int(round(reynolds_lattice * viscosity_lattice / velocity_lattice))     # lattices in x-direction
M = int(round(N / aspect))                                                  # lattices in y-direction (want dx = dy)

# weights for D2Q9
w = np.zeros(9)
w[0] = 4./9.
w[1:5] = 1./9.
w[5:9] = 1./36.

# grids MxNx9
f = np.zeros((9,M,N))       # distribution function
#fnext = np.zeros((9,M,N))

rho = np.zeros((M,N))       # density
u = np.zeros((2,M,N))       # velocity

# boundary condition
u0 = np.array((velocity_lattice,0.))
u[0,0,1:-1] = velocity_lattice     # load BC into north boundary

# initial condition
rho[:,:] = rho0

# unit vectors for D2Q9
c = np.array(((0,0),(1,0),(0,1),(-1,0),(0,-1),(1,1),(-1,1),(-1,-1),(1,-1)))

# sim length
nsteps = 40000


# for statistics collection
rhosum_bound = np.zeros((5,nsteps+1))     # 0-4 Tot E N W S
rhosum_sect = np.zeros((2,nsteps+1))      # 0-1 horz vert

rhosum_bound[0,0] = rho.sum()           # total
rhosum_bound[1,0] = rho[:,-1].sum()     # E
rhosum_bound[2,0] = rho[0,:].sum()      # N
rhosum_bound[3,0] = rho[:,0].sum()      # W
rhosum_bound[4,0] = rho[-1,:].sum()     # S

rhosum_sect[0,0] = rho[M/2,:].sum()   # horz
rhosum_sect[1,0] = rho[:,N/2].sum()   # vert


print 'Reynolds:'+str(reynolds_lattice)
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


