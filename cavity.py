from __future__ import division

__author__ = 'mh'

# D2Q9 Lattice Boltzmann

import numpy as np

# collision

def equilibrium(rho,u,w,c):
    # rho[M,N]
    # u[2,M,N]  first dimension is u and v
    # c[9,2]    first dimension is k, second is cx and cy
    m = rho.shape[0]
    n = rho.shape[1]

    # cu = 3.0 * np.dot(c,u.transpose(1,0,2))
    # usqr = 3./2.*(u[0]**2+u[1]**2)
    # feq = zeros((q,nx,ny))
    # for i in range(q): feq[i,:,:] = rho*t[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
    # return feq

    cTu = np.zeros((9,m,n))
    for i in range(9):
        cTu[i] = c[i,0]*u[0] + c[i,1]*u[1]

    uTu = u[0]**2 + u[1]**2

    feq = np.zeros((9,m,n))

    for k in range(0,9):
        cx = c[k,0]
        cy = c[k,1]
        wk = w[k]

        #ckTuij = 3.*((cx * u[0]) + (cy * u[1]))
        #uijTuij = u[0]**2 + u[1]**2

        # feq[k,:,:] = (rho * wk * (1. + (3. * ckTuij) + (4.5 * np.square(ckTuij)) - (1.5 * uijTuij) ))
        feq[k,:,:] = (rho * wk * (1. + (3. * cTu[k]) + (0.5 * (cTu[k] ** 2)) - (1.5 * uTu) ))

    return feq

def otherEquilibrium(rho,u,w,c):

    nx = rho.shape[0]
    ny = rho.shape[1]

    cu = 3.0 * np.dot(c,u.transpose(1,0,2))
    #cu = np.dot(c,u.transpose(1,0,2))

    CU = np.zeros((9,nx,ny))
    for i in range(9):
        CU[i] = 3.*(c[i,0]*u[0] + c[i,1]*u[1])

    print 'diff cu-CU:'
    print (cu-CU).sum()

    usqr = 3./2.*(u[0]**2+u[1]**2)
    feq = np.zeros((9,nx,ny))

    for i in range(9):
        feq[i,:,:] = rho*w[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
    return feq

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
        uijTuij = u[0]**2 + u[1]**2

        # update fk with feqk
        #f[k,:,:] *= (1. - omega)
        f[k,:,:] = f[k,:,:]*(1.-omega) + (rho * wk * (1. + (3. * ckTuij) + (4.5 * np.square(ckTuij)) - (1.5 * uijTuij))) * omega

        return f

# streaming

def streaming(f):
    # k = f.shape[0]
    # m = f.shape[1]
    # n = f.shape[2]

    # f[1,:,1:] = f[1,:,:(n-1)]  # stream one step E
    # f[2,:(m-1),:] = f[2,1:,:]  # stream one step N
    # f[3,:,:(n-1)] = f[3,:,1:]  # stream one step W
    # f[4,1:,:] = f[4,:(m-1),:]  # stream one step S
    #
    # f[5,:(m-1),1:] = f[5,1:,:(n-1)]     # stream one step NE
    # f[6,:(m-1),:(n-1)] = f[6,1:,1:]     # stream one step NW
    # f[7,1:,:(n-1)] = f[7,:(m-1),1:]     # stream one step SW
    # f[8,1:,1:] = f[8,:(m-1),:(n-1)]     # stream one step SE

    X = 1
    Y = 0

    E = 1
    N = -1
    W = -1
    S = 1

    f[1,:,:] = np.roll(f[1,:,:],E,axis=X)  # E
    f[2,:,:] = np.roll(f[2,:,:],N,axis=Y)  # N
    f[3,:,:] = np.roll(f[3,:,:],W,axis=X)  # W
    f[4,:,:] = np.roll(f[4,:,:],S,axis=Y)  # S

    f[5,:,:] = np.roll(np.roll(f[5,:,:],E,axis=X),N,axis=Y)  # NE
    f[6,:,:] = np.roll(np.roll(f[6,:,:],N,axis=Y),W,axis=X)  # NW
    f[7,:,:] = np.roll(np.roll(f[7,:,:],W,axis=X),S,axis=Y)  # SW
    f[8,:,:] = np.roll(np.roll(f[8,:,:],S,axis=Y),E,axis=X)  # SE

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

    rhoN = f[0,0,:] + f[1,0,:] + f[3,0,:] + (2 * (f[2,0,:] + f[5,0,:] + f[6,0,:]))

    s = (1./6) * rhoN * u0[0]
    #s = 0.5 * rhoN * u0[0]    # TODO is this correct? or should it be 1/6?

    f[4,0,:] = f[2,0,:]
    #f[7,0,:] = f[5,0,:] + 0.5 * (f[1,0,:] - f[3,0,:]) - s # TODO is it correct to assume (f1-f3)=0 ?
    #f[8,0,:] = f[6,0,:] + 0.5 * (f[3,0,:] - f[1,0,:]) + s
    f[7,0,:] = f[5,0,:] - s
    f[8,0,:] = f[6,0,:] + s

    return f


# calculate density and velocity component

def density(f):

     # f[9,M,N]  first dimension is k
    #rho = np.sum(f,axis=0)    # new rho
    rho = f[0,:,:]+f[1,:,:]+f[2,:,:]+f[3,:,:]+f[4,:,:]+f[5,:,:]+f[6,:,:]+f[7,:,:]+f[8,:,:]
    rho[0,:] = f[0,0,:] + f[1,0,:] + f[3,0,:] + (2.0 * (f[2,0,:] + f[5,0,:] + f[6,0,:]))

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

    #u /= rho

    return u/rho



