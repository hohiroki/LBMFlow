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

    cTu = np.zeros((9,m,n))
    for k in range(9):
        cTu[k] = 3. * (c[k,0]*u[0] + c[k,1]*u[1])

    uTu = 1.5 * (u[0]**2 + u[1]**2)

    feq = np.zeros((9,m,n))

    for k in range(0,9):

        feq[k] = rho * w[k] * (1. + cTu[k] + (0.5 * np.square(cTu[k])) - uTu )

    return feq

def otherEquilibrium(rho,u,w,c):

    nx = rho.shape[0]
    ny = rho.shape[1]

    cu = 3.0 * np.dot(c,u.transpose(1,0,2))
    #cu = np.dot(c,u.transpose(1,0,2))

    usqr = 3./2.*(u[0]**2+u[1]**2)
    feq = np.zeros((9,nx,ny))

    for i in range(9):
        feq[i,:,:] = rho * w[i] * (1.+cu[i]+0.5*cu[i]**2-usqr)

    return feq

def otherEquilibrium2(rho,u,w,c):

    nx = rho.shape[0]
    ny = rho.shape[1]

    cu = 3.0 * np.dot(c,u.transpose(1,0,2))
    #cu = np.dot(c,u.transpose(1,0,2))

    #usqr = 3./2.*(u[0]**2+u[1]**2)
    #usqr = 3./2.*( np.dot(u[0],u[0]) + np.dot(u[1],u[1]) )
    usqr = np.square(u).sum(0)
    feq = np.zeros((9,nx,ny))

    for i in range(9):
        feq[i,:,:] = rho * w[i] * (1.+cu[i]+0.5*cu[i]**2-1.5*usqr)

    return feq

def otherEquilibrium3(rho,u,w,c):

    nx = rho.shape[0]
    ny = rho.shape[1]

    tcu = 3.0 * np.dot(c,u.transpose(1,0,2))
    cu = np.dot(c,u.transpose(1,0,2))
    #cusqr = 4.5*(np.dot(c,u.transpose(1,0,2))**2)
    cusqr = np.dot(cu,cu)

    usqr = 1.5*(u[0]**2+u[1]**2)
    feq = np.zeros((9,nx,ny))

    for i in range(9):
        feq[i,:,:] = rho * w[i] * (1. + cu[i] + cusqr[i] - usqr)

    return feq

# def collisionTwo(f,feq,omega):
#
#     # f[9,M,N]  first dimension is k
#     # u[2,M,N]  first dimension is u and v
#     # c[9,2]    first dimension is k, second is cx and cy
#     # rho[M,N]
#     # w[9]
#     # omega[scalar]
#
#     # calculate one k at a time, use matrices
#     m = f.shape[1]
#     n = f.shape[2]
#
#     newf = np.zeros((9,m,n))
#
#     oneminomega = 1. - omega
#
#     for k in range(9):
#
#         # # helpers for feqk
#         # cx = c[k,0]
#         # cy = c[k,1]
#         # wk = w[k]
#         #
#         # ckTuij = (cx * u[0,:,:]) + (cy * u[1,:,:])
#         # uijTuij = u[0]**2 + u[1]**2
#
#         # update fk with feqk
#         #f[k,:,:] *= (1. - omega)
#         #f[k,:,:] = f[k,:,:]*(1.-omega) + (rho * wk * (1. + (3. * ckTuij) + (4.5 * np.square(ckTuij)) - (1.5 * uijTuij))) * omega
#
#         newf[k] = (f[k]*oneminomega) + (feq[k] * omega)
#
#     return newf

def collision(f,feq,omega):
    oneminomega = 1. - omega
    return (f*oneminomega) + (feq*omega)

# def otherCollision(f,feq,omega):
#     return f - omega * (f - feq)

# streaming

# def streamingThree(f):
#
#     m = f.shape[1]
#     n = f.shape[2]
#
#     newf = np.zeros((9,m,n))
#
#
#     k = f.shape[0]
#     m = f.shape[1]
#     n = f.shape[2]
#
#     newf[0] = f[0]
#
#     newf[1,:,1:] = f[1,:,:(n-1)]  # stream one step E
#     newf[2,:(m-1),:] = f[2,1:,:]  # stream one step N
#     newf[3,:,:(n-1)] = f[3,:,1:]  # stream one step W
#     newf[4,1:,:] = f[4,:(m-1),:]  # stream one step S
#
#     newf[5,:(m-1),1:] = f[5,1:,:(n-1)]     # stream one step NE
#     newf[6,:(m-1),:(n-1)] = f[6,1:,1:]     # stream one step NW
#     newf[7,1:,:(n-1)] = f[7,:(m-1),1:]     # stream one step SW
#     newf[8,1:,1:] = f[8,:(m-1),:(n-1)]     # stream one step SE
#
#     return newf

def streaming(f):

    # f[0] = f[0]

    f[1,:,1:] = f[1,:,:-1]      # stream one step E
    f[2,:-1,:] = f[2,1:,:]      # stream one step N
    f[3,:,:-1] = f[3,:,1:]      # stream one step W
    f[4,1:,:] = f[4,:-1,:]      # stream one step S

    f[5,:-1,1:] = f[5,1:,:-1]           # stream one step NE
    f[6,:-1,:-1] = f[6,1:,1:]           # stream one step NW
    f[7,1:,:-1] = f[7,:-1,1:]           # stream one step SW
    f[8,1:,1:] = f[8,:-1,:-1]           # stream one step SE

    return f



# def streamingTwo(f):
#
#     m = f.shape[1]
#     n = f.shape[2]
#
#     newf = np.zeros((9,m,n))
#
#     X = 1
#     Y = 0
#
#     E = 1
#     N = -1
#     W = -1
#     S = 1
#
#     newf[0] = f[0]
#
#     newf[1,:,:] = np.roll(f[1,:,:],E,axis=X)  # E
#     newf[2,:,:] = np.roll(f[2,:,:],N,axis=Y)  # N
#     newf[3,:,:] = np.roll(f[3,:,:],W,axis=X)  # W
#     newf[4,:,:] = np.roll(f[4,:,:],S,axis=Y)  # S
#
#     newf[5,:,:] = np.roll(np.roll(f[5,:,:],E,axis=X),N,axis=Y)  # NE
#     newf[6,:,:] = np.roll(np.roll(f[6,:,:],N,axis=Y),W,axis=X)  # NW
#     newf[7,:,:] = np.roll(np.roll(f[7,:,:],W,axis=X),S,axis=Y)  # SW
#     newf[8,:,:] = np.roll(np.roll(f[8,:,:],S,axis=Y),E,axis=X)  # SE
#
#     return newf

# calculate distribution function at boundaries

def westBoundary(f):
    # f[9,M,N]  first dimension is k

    # bounce back
    f[1,:,0] = f[3,:,0]
    f[5,:,0] = f[7,:,0]
    f[8,:,0] = f[6,:,0]

    return f

def eastBoundary(f):
    # f[9,M,N]  first dimension is k

    # bounce back
    f[3,:,-1] = f[1,:,-1]
    f[7,:,-1] = f[5,:,-1]
    f[6,:,-1] = f[8,:,-1]

    return f

def southBoundary(f):
    # f[9,M,N]  first dimension is k

    # bounce back
    f[2,-1,:] = f[4,-1,:]
    f[5,-1,:] = f[7,-1,:]
    f[6,-1,:] = f[8,-1,:]

    return f

def northBoundary(f,u0):

    rhoN = f[0,0,1:-1] + f[1,0,1:-1] + f[3,0,1:-1] + (2 * (f[2,0,1:-1] + f[5,0,1:-1] + f[6,0,1:-1]))

    s = (1./6) * rhoN * u0[0]   # we call this alt A
    #s = 0.5 * rhoN * u0[0]    # alt B

    f[4,0,1:-1] = f[2,0,1:-1]
    #f[7,0,1:-1] = f[5,0,1:-1] + 0.5 * (f[1,0,1:-1] - f[3,0,1:-1]) - s # TODO is it correct to assume (f1-f3)=0 ?
    #f[8,0,1:-1] = f[6,0,1:-1] + 0.5 * (f[3,0,1:-1] - f[1,0,1:-1]) + s
    #f[7,0,1:-1] = f[5,0,1:-1] - s
    #f[8,0,1:-1] = f[6,0,1:-1] + s

    f[7,0,1:-1] = f[5,0,1:-1] + (1./6) * (f[1,0,1:-1] - f[3,0,1:-1]) - s # altC
    f[8,0,1:-1] = f[6,0,1:-1] + (1./6) * (f[3,0,1:-1] - f[1,0,1:-1]) + s

    return f,rhoN


# calculate density and velocity component

def density(f,rhoN):

    rho = np.sum(f,axis=0)
    rhoNN = f[0,0,1:] + f[1,0,1:] + f[3,0,1:] + (2 * (f[2,0,1:] + f[5,0,1:] + f[6,0,1:]))
    rho[0,1:]=rhoNN

    #rho = np.sum(f,axis=0)
    #rho[0,1:-1]=rhoN

    return rho

def fasterVelocity(f,rho,c):

    # vel = np.dot(c.transpose(), f.transpose((1,0,2)))/rho
    # vel[:,0,1:] = u[:,0,1:]
    # vel[:,:,0] = u[:,:,0]

    return np.dot(c.transpose(), f.transpose((1,0,2)))/rho


def velocity(f,rho,c):

    m = rho.shape[0]
    n = rho.shape[1]

    fcx = np.zeros((9,m,n))
    fcy = np.zeros((9,m,n))
    for k in range(9):
        # fcx[k] = f[k]*c[k,0]
        # fcy[k] = f[k]*c[k,1]
        fcx[k] = np.dot(f[k],c[k,0])
        fcy[k] = np.dot(f[k],c[k,1])

    return np.array((fcx.sum(0),fcy.sum(0)))/rho


# def streamfunction(u,rho):
#
#     m = rho.shape[0]
#     n = rho.shape[1]
#     v = u[1]
#
#     strf = np.zeros((m,n))
#
#     for i in range(n):
#         rhoav = 0.5*(rho[-1,i])
#
#         if i != 0:
#             rhoav = 0.5*(rho[-1,i-1]+rho[-1,i])
#             strf[-1,i] = strf[-1,i-1]-rhoav*0.5*(v[-1,i-1]+v[-1,i])
#         for j in range(1,m):
#
#
#     rhoav = rho[-1,:].copy()
#     rhoav[1:] += rho[-1,:-1]
#     rhoav *= 0.5
#
#     s = 0.5*rhoav*(v[-1,:-1] + v[-1,1:])
#
#     strf[-1:,1:]
#
#        t = strf[-1,1:].copy()
#
#
#
#
#     return strf