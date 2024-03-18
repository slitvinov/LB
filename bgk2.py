import numpy as np

# 0 = rest particles, 1-4 = velocity 1, 4-9 = velocity sqrt(2)

def equil():
    for j in range(ny + 2):
        for i in range(nx + 2):
            rl = rho[i, j]
            ul = u[i, j] / cs2
            vl = v[i, j] / cs2
            uv = ul * vl
            usq = u[i, j] * u[i, j]
            vsq = v[i, j] * v[i, j]
            sumsq = (usq + vsq) / cs22
            sumsq2 = sumsq * (1.0 - cs2) / cs2
            u2 = usq / cssq
            v2 = vsq / cssq
            feq[0, i, j] = w0 * (1.0 - sumsq)

            feq[1, i, j] = w1 * (1.0 - sumsq + u2 + ul)
            feq[2, i, j] = w1 * (1.0 - sumsq + v2 + vl)
            feq[3, i, j] = w1 * (1.0 - sumsq + u2 - ul)
            feq[4, i, j] = w1 * (1.0 - sumsq + v2 - vl)

            feq[5, i, j] = w2 * (1.0 + sumsq2 + ul + vl + uv)
            feq[6, i, j] = w2 * (1.0 + sumsq2 - ul + vl - uv)
            feq[7, i, j] = w2 * (1.0 + sumsq2 - ul - vl + uv)
            feq[8, i, j] = w2 * (1.0 + sumsq2 + ul - vl - uv)


nx = 128
ny = 64
npop = 9
nsteps = 101
omega = 1.5
iforce = True
rho0 = 1
u0 = 0.1
v0 = 0
uf = 0.1
iobst = True
nobst = 8

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

u = np.full((nx + 2, ny + 2), u0, np.float64)
v = np.full((nx + 2, ny + 2), v0, np.float64)
rho = np.full((nx + 2, ny + 2), rho0, np.float64)
feq = np.empty((npop, nx + 2, ny + 2), np.float64)
equil()
f = np.copy(feq)
for istep in range(1, nsteps + 1):
    # inlet
    f[1, 0, :] = f[1, nx, :]
    f[5, 0, :] = f[5, nx, :]
    f[8, 0, :] = f[8, nx, :]

    # outlet
    f[3, nx + 1, :] = f[3, 1, :]
    f[6, nx + 1, :] = f[6, 1, :]
    f[7, nx + 1, :] = f[7, 1, :]

    # bounce back (north)
    for i in range(1, nx + 1):
        f[4, i, ny + 1] = f[2, i, ny]
        f[8, i, ny + 1] = f[6, i + 1, ny]
        f[7, i, ny + 1] = f[5, i - 1, ny]

    # bounce back (source)
    for i in range(1, nx + 1):
        f[2, i, 0] = f[4, i, 1]
        f[6, i, 0] = f[8, i - 1, 1]
        f[5, i, 0] = f[7, i + 1, 1]

    # corners
    f[8, 0, ny + 1] = f[6, 1, ny]
    f[5, 0, 0] = f[7, 1, 1]
    f[7, nx + 1, ny + 1] = f[5, nx, ny]
    f[6, nx + 1, 0] = f[8, nx, 1]

    # move
    for j in range(ny, 0, -1):
        for i in range(1, nx + 1):
            f[2, i, j] = f[2, i, j - 1]
            f[6, i, j] = f[6, i + 1, j - 1]
    for j in range(ny, 0, -1):
        for i in range(nx, 0, -1):
            f[1, i, j] = f[1, i - 1, j]
            f[5, i, j] = f[5, i - 1, j - 1]

    for j in range(1, ny + 1):
        for i in range(nx, 0, -1):
            f[4, i, j] = f[4, i, j + 1]
            f[8, i, j] = f[8, i - 1, j + 1]

    for j in range(1, ny + 1):
        for i in range(1, nx + 1):
            f[3, i, j] = f[3, i + 1, j]
            f[7, i, j] = f[7, i + 1, j + 1]

    # hydrovar
    for j in range(1, ny + 1):
        for i in range(1, nx + 1):
            rho[i, j] = f[1, i, j] + f[2, i, j] + f[3, i, j] + f[4, i, j] + f[
                5, i, j] + f[6, i, j] + f[7, i, j] + f[8, i, j] + f[0, i, j]
            rhoi = 1. / rho[i, j]
            u[i, j] = (f[1, i, j] - f[3, i, j] + f[5, i, j] - f[6, i, j] -
                       f[7, i, j] + f[8, i, j]) * rhoi
            v[i, j] = (f[5, i, j] + f[2, i, j] + f[6, i, j] - f[7, i, j] -
                       f[4, i, j] - f[8, i, j]) * rhoi
    equil()
    # collision step
    f = f * (1.0 - omega) + omega * feq
    
    if iforce:
        f[1, ::] += fpois
        f[5, ::] += fpois
        f[8, ::] += fpois
        
        f[3, ::] -= fpois
        f[6, ::] -= fpois
        f[7, ::] -= fpois
        
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
            file.write(u[1:nx + 1, 1:ny + 1].tobytes("F"))
            file.write(v[1:nx + 1, 1:ny + 1].tobytes("F"))
            file.write(rho[1:nx + 1, 1:ny + 1].tobytes("F"))
            vort = np.roll(u, [0, 1]) - np.roll(u, [0, -1]) - np.roll(
                v, [1, 0]) + np.roll(v, [-1, 0])
            file.write(vort[1:nx + 1, 1:ny + 1].tobytes("F"))
print(np.var(u), np.var(v), np.var(rho))
