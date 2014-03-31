from __future__ import division

__author__ = 'mh'

# D2Q9 Lattice Boltzmann

import numpy as np

# input data

# fluid physics
viscosity = 1.2*1e-3    # kinematic viscosity m2/s (oil at 15 deg C)

velocity = 6.0           # moving lid m/s

height = 0.2     # height of cavity m
width = 0.2      # width of cavity m

reynolds = velocity * height / viscosity    # macroscale reynolds number

# LBM values
aspect = width / height         # aspect ratio of domain

reynolds_lattice = reynolds     # keep same reynolds
viscosity_lattice = 0.01        # choose so we keep same reynolds
velocity_lattice = 0.1          # choose so we keep same reynolds (here we ignore the actual velocity)

omega = 1./(3. * viscosity_lattice + 0.5)
#omega = (3. * velocity_lattice) + 0.5

N = int(round(reynolds_lattice * viscosity_lattice / velocity_lattice))     # lattices in x-direction
M = int(round(N / aspect))                                                  # lattices in y-direction (want dx = dy)

# weights for D2Q9
w = np.zeros(9)
w[0] = 4./9.
w[1:5] = 1./9.
w[5:9] = 1./36.

# grids MxNx9
f = np.zeros((9,M,N))       # distribution function
rho = np.zeros((M,N))       # density
u = np.zeros((2,M,N))       # velocity

# boundary condition
u0 = np.array((velocity_lattice,0.))
#u[0,0,1:-1] = velocity_lattice     # load BC into north boundary

# initial condition
rho[:,:] = 5.

# unit vectors for D2Q9
c = np.array(((0,0),(1,0),(0,1),(-1,0),(0,-1),(1,1),(-1,1),(-1,-1),(1,-1)))




