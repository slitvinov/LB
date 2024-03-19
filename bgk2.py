#import numpy as np
import torch as np
import copy


# 0 = rest particles, 1-4 = velocity 1, 4-9 = velocity sqrt(2)
def equil():
    ul = u / cs2
    vl = v / cs2
    uv = ul * vl
    usq = u * u
    vsq = v * v
    sumsq = (usq + vsq) / cs22
    sumsq2 = sumsq * (1.0 - cs2) / cs2
    u2 = usq / cssq
    v2 = vsq / cssq
    feq[0] = w0 * (1.0 - sumsq)

    feq[1] = w1 * (1.0 - sumsq + u2 + ul)
    feq[2] = w1 * (1.0 - sumsq + v2 + vl)
    feq[3] = w1 * (1.0 - sumsq + u2 - ul)
    feq[4] = w1 * (1.0 - sumsq + v2 - vl)

    feq[5] = w2 * (1.0 + sumsq2 + ul + vl + uv)
    feq[6] = w2 * (1.0 + sumsq2 - ul + vl - uv)
    feq[7] = w2 * (1.0 + sumsq2 - ul - vl + uv)
    feq[8] = w2 * (1.0 + sumsq2 + ul - vl - uv)


nx = 128 * 8
ny = 64 * 8
nsteps = 10001
omega = 1.5
iforce = True
rho0 = 1
u0 = 0.1
v0 = 0
uf = 0.1
iobst = True
nobst = 8 * 8

# lattice weights
w0 = 4.0 / 9.0
w1 = 1.0 / 9.0
w2 = 1.0 / 36.0

# sound - speed and related constants
cs2 = 1.0 / 3.0
cs22 = 2.0 * cs2
cssq = 2.0 / 9.0
visc = (1.0 / omega - 0.5) * cs2
rey = u0 * ny / visc

# applied force (based on Stokes problem)
fpois = 8.0 * visc * uf / ny / ny
fpois = rho0 * fpois / 6.

u = np.full((nx + 2, ny + 2), u0, dtype=np.float64)
v = np.full((nx + 2, ny + 2), v0, dtype=np.float64)
rho = np.full((nx + 2, ny + 2), rho0, dtype=np.float64)
feq = np.empty((9, nx + 2, ny + 2), dtype=np.float64)
equil()
f = feq[:]
for istep in range(1, nsteps + 1):
    # inlet
    f[1, 0, 1:ny + 1] = f[1, nx, 1:ny + 1]
    f[5, 0, 1:ny + 1] = f[5, nx, 1:ny + 1]
    f[8, 0, 1:ny + 1] = f[8, nx, 1:ny + 1]

    # outlet
    f[3, nx + 1, 1:ny + 1] = f[3, 1, 1:ny + 1]
    f[6, nx + 1, 1:ny + 1] = f[6, 1, 1:ny + 1]
    f[7, nx + 1, 1:ny + 1] = f[7, 1, 1:ny + 1]

    # bounce back (north)
    f[4, 1:nx + 1, ny + 1] = f[2, 1:nx + 1, ny]
    f[8, 1:nx + 1, ny + 1] = f[6, 2:, ny]
    f[7, 1:nx + 1, ny + 1] = f[5, :nx, ny]

    # bounce back (source)
    f[2, 1:nx + 1, 0] = f[4, 1:nx + 1, 1]
    f[6, 1:nx + 1, 0] = f[8, :nx, 1]
    f[5, 1:nx + 1, 0] = f[7, 2:, 1]

    # corners
    f[8, 0, ny + 1] = f[6, 1, ny]
    f[5, 0, 0] = f[7, 1, 1]
    f[7, nx + 1, ny + 1] = f[5, nx, ny]
    f[6, nx + 1, 0] = f[8, nx, 1]

    # move
    f[2, 1:nx + 1, 1:ny + 1] = copy.deepcopy(f[2, 1:nx + 1, :ny])
    f[6, 1:nx + 1, 1:ny + 1] = copy.deepcopy(f[6, 2::, :ny])
    f[1, 1:nx + 1, 1:ny + 1] = copy.deepcopy(f[1, :nx, 1:ny + 1])
    f[5, 1:nx + 1, 1:ny + 1] = copy.deepcopy(f[5, :nx, :ny])
    f[4, 1:nx + 1, 1:ny + 1] = copy.deepcopy(f[4, 1:nx + 1, 2:])
    f[8, 1:nx + 1, 1:ny + 1] = copy.deepcopy(f[8, :nx, 2:])
    f[3, 1:nx, 1:ny + 1] = copy.deepcopy(f[3, 2:nx + 1, 1:ny + 1])
    f[7, 1:nx, 1:ny + 1] = copy.deepcopy(f[7, 2:nx + 1, 2:])

    # hydrovar
    rho[:] = f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[0]
    u = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) / rho
    v = (f[5] + f[2] + f[6] - f[7] - f[4] - f[8]) / rho
    equil()
    # collision step
    f = f * (1.0 - omega) + omega * feq

    if iforce:
        f[1] += fpois
        f[5] += fpois
        f[8] += fpois

        f[3] -= fpois
        f[6] -= fpois
        f[7] -= fpois

    if iobst:
        i = nx // 4
        jbot = ny // 2 - nobst // 2
        jtop = ny // 2 + nobst // 2 + 1
        for j in range(ny // 2 - nobst // 2, ny // 2 + nobst // 2 + 2):
            f[1, i, j] = f[3, i + 1, j]
            f[5, i, j] = f[7, i + 1, j + 1]
            f[8, i, j] = f[6, i + 1, j - 1]
            f[3, i, j] = f[1, i - 1, j]
            f[7, i, j] = f[5, i - 1, j]
            f[6, i, j] = f[8, i - 1, j]
        f[2, i, jtop] = f[4, i, jtop + 1]
        f[6, i, jtop] = f[8, i - 1, jtop + 1]
        f[4, i, jbot] = f[2, i, jbot - 1]
        f[7, i, jbot] = f[5, i - 1, jbot - 1]
    if istep % 100 == 0:
        path = "bgk.%08d.raw" % istep
        with open(path, "wb") as file:
            vort = np.roll(u, [0, 1], [0, 1]) - np.roll(
                u, [0, -1], [0, 1]) - np.roll(v, [1, 0], [0, 1]) + np.roll(
                    v, [-1, 0], [0, 1])
            for field in u, v, rho, vort:
                if hasattr(field, "numpy"):
                    field = field.detach().numpy()
                file.write(field[1:nx + 1, 1:ny + 1].tobytes("F"))
        print(np.var(u[1:nx + 1, 1:ny + 1]), np.var(v[1:nx + 1, 1:ny + 1]),
              np.var(rho[1:nx + 1, 1:ny + 1]))
