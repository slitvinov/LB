import numpy as np

nx = 400
ny = 100
obst_x = nx // 5 + 1
obst_y = ny // 2 + 3
obst_r = ny // 10 + 1
uMax = 0.1
Re = 100
nu = uMax * 2. * obst_r / Re
omega = 1. / (3 * nu + 1. / 2.)
maxT = 400000
tPlot = 50
t = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36])
cx = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1])
cy = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1])
opp = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])
col = np.arange(1, ny - 1)
x, y = np.meshgrid(np.arange(ny), np.arange(nx))
obst = (x - obst_x) ** 2 + (y - obst_y) ** 2 <= obst_r ** 2
obst[:, [0, ny - 1]] = 1
bbRegion = np.where(obst)
L = ny - 2
y_phys = y - 1.5
u = 4 * uMax / (L * L) * (y_phys * L - y_phys * y_phys)
v = np.zeros((nx, ny))
rho = 1
f0 = np.zeros((9, nx, ny))

for i in range(9):
    cu = 3 * (cx[i] * u + cy[i] * v)
    f0[i, :, :] = rho * t[i] * (1 + cu + 1 / 2 * (cu * cu) - 3 / 2 * (u ** 2 + v ** 2))

for cycle in range(maxT):
    rho = np.sum(f0, axis=0)
    u = np.sum(cx[:, np.newaxis, np.newaxis] * f0, axis=0) / rho
    v = np.sum(cy[:, np.newaxis, np.newaxis] * f0, axis=0) / rho

    # MACROSCOPIC (DIRICHLET) BOUNDARY CONDITIONS
    # inlet: Poiseuille profile
    y_phys = col - 1.5
    u[0, col] = 4 * uMax / (L * L) * (y_phys * L - y_phys * y_phys)
    v[0, col] = 0
    rho[0, col] = 1 / (1 - u[0, col]) * (np.sum(f0[[0, 2, 4], 0][:, col], axis=0)
                                         + 2 * np.sum(f0[[3, 6, 7], 0][:, col], axis=0))

    # outlet: Constant pressure
    rho[nx - 1, col] = 1
    u[nx - 1, col] = -1 + 1 / (rho[nx - 1, col]) * (np.sum(f0[[0, 2, 4], nx - 1][:, col], axis=0)
                                                         + 2 * np.sum(f0[[1, 5, 8], nx - 1][:, col], axis=0))
    v[nx - 1, col] = 0

    # MICROSCOPIC BOUNDARY CONDITIONS: inlet (Zou/He BC)
    f0[1, 0, col] = f0[3, 0, col] + 2 / 3 * rho[0, col] * u[0, col]
    f0[5, 0, col] = f0[7, 0, col] + 1 / 2 * (f0[4, 0, col] - f0[2, 0, col]) \
                    + 1 / 2 * rho[0, col] * v[0, col] \
                    + 1 / 6 * rho[0, col] * u[0, col]
    f0[8, 0, col] = f0[6, 0, col] + 1 / 2 * (f0[2, 0, col] - f0[4, 0, col]) \
                    - 1 / 2 * rho[0, col] * v[0, col] \
                    + 1 / 6 * rho[0, col] * u[0, col]

    # MICROSCOPIC BOUNDARY CONDITIONS: outlet (Zou/He BC)
    f0[3, nx - 1, col] = f0[1, nx - 1, col] - 2 / 3 * rho[nx - 1, col] * u[nx - 1, col]
    f0[7, nx - 1, col] = f0[5, nx - 1, col] + 1 / 2 * (f0[2, nx - 1, col] - f0[4, nx - 1, col]) \
                         - 1 / 2 * rho[nx - 1, col] * v[nx - 1, col] \
                         - 1 / 6 * rho[nx - 1, col] * u[nx - 1, col]
    f0[6, nx - 1, col] = f0[8, nx - 1, col] + 1 / 2 * (f0[4, nx - 1, col] - f0[2, nx - 1, col]) \
                         + 1 / 2 * rho[nx - 1, col] * v[nx - 1, col] \
                         - 1 / 6 * rho[nx - 1, col] * u[nx - 1, col]

    # Calculate equilibrium distribution
    feq = np.zeros_like(f0)
    for i in range(9):
        cu = 3 * (cx[i] * u + cy[i] * v)
        feq[i, :, :] = rho * t[i] * (1 + cu + 1 / 2 * (cu * cu) - 3 / 2 * (u ** 2 + v ** 2))

    # Collision step
    f1 = f0 - omega * (f0 - feq)

    # Bounce-back boundary condition
    for i in range(9):
        f1[i, bbRegion[0], bbRegion[1]] = f0[opp[i], bbRegion[0], bbRegion[1]]

    # Streaming step
    for i in range(9):
        f0[i, :, :] = np.roll(f1[i, :, :], (cx[i], cy[i]), axis=(0, 1))

    if cycle % tPlot == 0:
        path = "cyl.%09d.raw" % cycle
        with open(path, "wb") as fid:
            for field in [u, v, rho, rho]:
                fid.write(field.tobytes())
        print(np.var(u), np.var(v), np.var(rho))
