import numpy as np
import matplotlib.pylab as plt

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
t = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36]
cx = [0, 1, 0, -1, 0, 1, -1, -1, 1]
cy = [0, 0, 1, 0, -1, 1, 1, -1, -1]
opp = [0, 3, 4, 1, 2, 7, 8, 5, 6]
col = np.arange(1, ny - 1)
y, x = np.meshgrid(np.arange(ny), np.arange(nx))
obst = (x - obst_x)**2 + (y - obst_y)**2 <= obst_r**2
obst[:, 0] = 1
obst[:, ny - 1] = 1
bb = np.where(obst)
L = ny - 2
y_phys = y - 1.5
u = 4 * uMax / (L * L) * (y_phys * L - y_phys * y_phys)
v = np.empty((nx, ny))
rho = np.empty((nx, ny))
rho0 = 1
f = np.empty((9, nx, ny))

for i in range(9):
    cu = 3 * (cx[i] * u + cy[i] * v)
    f[i, :, :] = rho0 * t[i] * (1 + cu + 1 / 2 * (cu * cu) - 3 / 2 *
                                (u**2 + v**2))

for cycle in range(maxT):
    rho[:] = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8]
    u[:] = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) / rho
    v[:] = (f[5] + f[2] + f[6] - f[7] - f[4] - f[8]) / rho

    # inlet: Poiseuille profile
    y_phys = col - 1.5
    u[0, col] = 4 * uMax / (L * L) * (y_phys * L - y_phys * y_phys)
    v[0, col] = 0
    rho[0, col] = 1 / (
        1 - u[0, col]) * (f[0][0, col] + f[2][0, col] + f[4][0, col] + 2 *
                          (f[3][0, col] + f[6][0, col] + f[7][0, col]))

    # outlet: Constant pressure
    rho[nx - 1, col] = 1
    u[nx - 1, col] = -1 + 1 / (rho[nx - 1, col]) * (
        f[0][nx - 1, col] + f[2][nx - 1, col] + f[4][nx - 1, col] + 2 *
        (f[1][nx - 1, col] + f[5][nx - 1, col] + f[8][nx - 1, col]))
    v[nx - 1, col] = 0

    # inlet (Zou/He BC)
    f[1][0, col] = f[3][0, col] + 2 / 3 * rho[0, col] * u[0, col]
    f[5][0, col] = f[7][0, col] + 1 / 2 * (f[4][0, col] - f[2][0, col]) \
                    + 1 / 2 * rho[0, col] * v[0, col] \
                    + 1 / 6 * rho[0, col] * u[0, col]
    f[8][0, col] = f[6][0, col] + 1 / 2 * (f[2][0, col] - f[4][0, col]) \
                    - 1 / 2 * rho[0, col] * v[0, col] \
                    + 1 / 6 * rho[0, col] * u[0, col]

    # outlet (Zou/He BC)
    f[3][nx - 1,
         col] = f[1][nx - 1, col] - 2 / 3 * rho[nx - 1, col] * u[nx - 1, col]
    f[7][nx - 1, col] = f[5][nx - 1, col] + 1 / 2 * (
        f[2][nx - 1, col] - f[4][nx - 1, col]) - 1 / 2 * rho[nx - 1, col] * v[
            nx - 1, col] - 1 / 6 * rho[nx - 1, col] * u[nx - 1, col]
    f[6][nx - 1, col] = f[8][nx - 1, col] + 1 / 2 * (
        f[4][nx - 1, col] - f[2][nx - 1, col]) + 1 / 2 * rho[nx - 1, col] * v[
            nx - 1, col] - 1 / 6 * rho[nx - 1, col] * u[nx - 1, col]

    feq = np.zeros_like(f)
    for i in range(9):
        cu = 3 * (cx[i] * u + cy[i] * v)
        feq[i][:] = rho * t[i] * (1 + cu + 1 / 2 * (cu * cu) - 3 / 2 *
                                  (u**2 + v**2))
    # Collision step
    f1 = f * (1.0 - omega) + omega * feq
    f1 = f - omega * (f - feq)

    # Bounce-back boundary condition
    for i in range(9):
        f1[i, bb[0], bb[1]] = f[opp[i], bb[0], bb[1]]

    # Streaming step
    for i in range(9):
        f[i, :, :] = np.roll(f1[i, :, :], (cx[i], cy[i]), axis=(0, 1))

    if cycle % tPlot == 0:
        path = "cyl.%09d.raw" % cycle
        with open(path, "wb") as fid:
            for field in [u, v, rho, rho]:
                fid.write(field.tobytes())
        print(np.var(u), np.var(v), np.var(rho))
