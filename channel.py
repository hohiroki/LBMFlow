__author__ = 'mh'

# input data


# physics
density = 1000.0        # density of water 1000.0 kg/m3
viscosity = 10**-6      # kinematic viscosity 10^-6 m2/s

vel_in = 0.02           # uniform velocity in on left side 0.2 m/s (characteristic vel macroscale)
channel_height = 0.02   # height of channel 2cm = 0.02m (characteristic length macroscale)
channel_length = 0.5    # length of channel 50cm = 0.5m

reynolds = vel_in * channel_height / viscosity    # macroscale reynolds number

# lattice values
aspect = channel_length / channel_height    # aspect ratio of domain

M = 40                  # lattices in y direction
N = int(aspect * M)     # lattices in x direction, we want deltax = deltay

vel_in_lattice = vel_in * 10    # choose so we can keep the other values same (sane?)
reynolds_lattice = reynolds     # keep this same as macroscale
viscosity_lattice = vel_in_lattice * M / reynolds_lattice   # use deltax = 1 so L = N (discretized channel_length)



# Main Loop


# calculate equilibrium distribution function feq


# collision


# streaming


# calculate distribution function at boundaries


# calculate density and velocity component


# End Main Loop


# output data
