from __future__ import division

__author__ = 'mh'

# D2Q9 Lattice Boltzmann

import numpy as np

# collision

def collision(f,rho,u,w,c,omega):

    # f[9,M,N]  first dimension is k
    # u[2,M,N]  first dimension is u and v
    # c[9,2]    first dimension is k, second is cx and cy
    # rho[M,N]
    # w[9]
    # omega[scalar]

    # calculate one k at a time, use matrices
    for k in range(0,9):

        # helpers for feqk
        cx = c[k,0]
        cy = c[k,1]
        wk = w[k]

        ckTuij = (cx * u[0,:,:]) + (cy * u[1,:,:])
        uijTuij = (np.square(u[0]) + np.square(u[1]))

        # update fk with feqk
        f[k,:,:] *= (1. - omega)
        f[k,:,:] += (rho * wk * (1. + 3. * ckTuij + 4.5 * np.square(ckTuij) - 1.5 * uijTuij)) * omega

        return f

# streaming

def streaming(f):
    # k = f.shape[0]
    m = f.shape[1]
    n = f.shape[2]

    f[1,:,1:] = f[1,:,:(n-1)]  # stream one step E
    f[2,:(m-1),:] = f[2,1:,:]  # stream one step N
    f[3,:,:(n-1)] = f[3,:,1:]  # stream one step W
    f[4,1:,:] = f[4,:(m-1),:]  # stream one step S

    f[5,:(m-1),1:] = f[5,1:,:(n-1)]     # stream one step NE
    f[6,:(m-1),:(n-1)] = f[6,1:,1:]     # stream one step NW
    f[7,1:,:(n-1)] = f[7,:(m-1),1:]     # stream one step SW
    f[8,1:,1:] = f[8,:(m-1),:(n-1)]     # stream one step SE

    return f


# calculate distribution function at boundaries

def westBoundary(f,u0):
    # f[9,M,N]  first dimension is k

    # bounce back
    f[1,:,0] = f[3,:,0]
    f[5,:,0] = f[7,:,0]
    f[8,:,0] = f[6,:,0]

    return f

def eastBoundary(f,u0):
    # f[9,M,N]  first dimension is k

    # bounce back
    f[3,:,-1] = f[1,:,-1]
    f[7,:,-1] = f[5,:,-1]
    f[6,:,-1] = f[8,:,-1]

    return f

def southBoundary(f,u0):
    # f[9,M,N]  first dimension is k

    # bounce back
    f[2,-1,:] = f[4,-1,:]
    f[5,-1,:] = f[7,-1,:]
    f[6,-1,:] = f[8,-1,:]

    return f

def northBoundary(f,u0):
    # f[9,M,N]  first dimension is k
    # u[2,M,N]  first dimension is u and v
    # u0[2]     BC x and y

    # ux0 = u[0,0,:]
    # f0 = f[0,0,:]
    # f1 = f[1,0,:]
    # f2 = f[2,0,:]
    # f3 = f[3,0,:]
    # f4 = f[4,0,:]
    # f5 = f[5,0,:]
    # f6 = f[6,0,:]
    # f7 = f[7,0,:]
    # f8 = f[8,0,:]

    # calculate three unknown fk and rhoN

    rhoN = f[0,0,:] + f[1,0,:] + f[3,0,:] + 2 * (f[2,0,:] + f[5,0,:] + f[6,0,:])

    #s = (1./6) * rhoN * u0[0]
    s = 0.5 * rhoN * u0[0]    # TODO is this correct? or should it be 1/6?

    f[4,0,:] = f[2,0,:]
    f[7,0,:] = f[5,0,:] + 0.5 * (f[1,0,:] - f[3,0,:]) - s # TODO is it correct to assume (f1-f3)=0 ?
    f[8,0,:] = f[6,0,:] + 0.5 * (f[3,0,:] - f[1,0,:]) + s
    # f[7,0,:] = f[5,0,:] - s
    # f[8,0,:] = f[6,0,:] + s

    return f


# calculate density and velocity component

def density(f,rho):

     # f[9,M,N]  first dimension is k
    rho = np.sum(f,axis=0)    # new rho

    return rho


def velocity(f,rho,u,c):

    # f[9,M,N]  first dimension is k
    # rho[M,N]
    # u[2,M,N]  first dimension is u and v
    # c[9,2]    first dimension is k, second is cx and cy

    u[:] = 0.   # set all velocities to 0

    # helpers for clarity
    # uij = u[0,:,:]
    # vij = u[1,:,:]

    cx = c[:,0]
    cy = c[:,1]

    for k in range(0,9):

        u[0,:,:] += f[k,:,:] * cx[k]
        u[1,:,:] += f[k,:,:] * cy[k]

    u /= rho

    return u



