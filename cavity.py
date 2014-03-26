__author__ = 'mh'

# input data


# fluid physics
viscosity = 1.2*1e-3    # kinematic viscosity m2/s (oil at 15 deg C)

velocity = 6.0           # moving lid m/s

height = 0.2     # height of cavity m
width = 0.2      # width of cavity m

reynolds = velocity * height / viscosity    # macroscale reynolds number

# LBM values
aspect = width / height    # aspect ratio of domain

reynolds_lattice = reynolds # keep same reynolds
viscosity_lattice = 0.01    # choose so we keep same reynolds
velocity_lattice = 0.1      # choose so we keep same reynolds

N = int(round(reynolds_lattice * viscosity_lattice / velocity_lattice)) # get lattices in x-direction
M = int(round(height * aspect)) # lattices in y-direction we want deltax = deltay



# Main Loop


# calculate equilibrium distribution function feq


# collision


# streaming


# calculate distribution function at boundaries


# calculate density and velocity component


# End Main Loop


# output data
