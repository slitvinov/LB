import numpy as np

nx = 128
ny = 64
npop = 9

nstep = 10001
nout = 2000
ndiag = 2000
omega = 1.5
iforce = True
rho0, u0, v0 = 1, 0.1, 0
uf = 0.1
ibobs = True
nobst = 8

#     lattice weights
w0 = 4.00 / 9.00
w1 = 1.0 / 9.00
w2 = 1.0 / 36.00

#     sound - speed and related constants
cs2 = 1.0 / 3.0
cs22 = 2.0 * cs2
cssq = 2.0 / 9.0

visc = (1.0 / omega - 0.5) * cs2
rey = u0 * ny / visc

#     Applied force(based on Stokes problem)
fpois = 8.0 * visc * uf / float(ny) / float(ny)

#     # of biased populations
fpois = rho0 * fpois / 6.

u = np.array((nx + 2, ny + 2))
v = np.array((nx + 2, ny + 2))
rho = np.array((nx + 2, ny + 2))
feq = np.array((npop, nx + 2, ny + 2))
f = np.array((npop, nx + 2, ny + 2))
