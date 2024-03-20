import numpy as np
# import torch as np

nx = 400
ny = 100
x0 = nx // 5 + 1
y0 = ny // 2 + 3
r0 = ny // 10 + 1
uMax = 0.1
Re = 100
nu = uMax * 2. * r0 / Re
omega = 1. / (3 * nu + 1. / 2.)
nsteps = 400000
tPlot = 500
t = 4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36
cx = 0, 1, 0, -1, 0, 1, -1, -1, 1
cy = 0, 0, 1, 0, -1, 1, 1, -1, -1
col = np.arange(1, ny - 1)
y, x = np.meshgrid(np.arange(ny), np.arange(nx))
obst = (x - x0)**2 + (y - y0)**2 <= r0**2
obst[:, 0] = 1
obst[:, ny - 1] = 1
xb, yb = np.where(obst)
xd, yd = np.where(~obst)
L = ny - 2
u = np.zeros((nx, ny))
v = np.zeros((nx, ny))
rho0 = 1
rho = np.full((nx, ny), rho0, dtype=np.float64)
f = [np.empty((nx, ny)) for i in range(9)]
feq = np.empty((nx, ny))

for i in range(9):
    c = 3 * (cx[i] * u + cy[i] * v)
    f[i][:] = rho0 * t[i] * (1 + c + 1 / 2 * (c * c) - 3 / 2 * (u**2 + v**2))

for cycle in range(nsteps):
    rho[:] = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8]
    u[:] = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) / rho
    v[:] = (f[5] + f[2] + f[6] - f[7] - f[4] - f[8]) / rho

    # inlet:
    y_phys = col - 0.5
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

    # collision and bounce-back boundary condition
    for i in range(9):
        c = 3 * (cx[i] * u + cy[i] * v)
        feq[:] = rho * t[i] * (1 + c + 1 / 2 * (c * c) - 3 / 2 * (u**2 + v**2))
        f[i][xd, yd] = f[i][xd, yd] * (1 - omega) + omega * feq[xd, yd]
    f[1][xb, yb], f[3][xb, yb] = f[3][xb, yb], f[1][xb, yb]
    f[2][xb, yb], f[4][xb, yb] = f[4][xb, yb], f[2][xb, yb]
    f[5][xb, yb], f[7][xb, yb] = f[7][xb, yb], f[5][xb, yb]
    f[6][xb, yb], f[8][xb, yb] = f[8][xb, yb], f[6][xb, yb]

    # Streaming step
    for i in range(9):
        f[i][:] = np.roll(f[i], (cx[i], cy[i]), axis=(0, 1))

    if cycle % tPlot == 0:
        print("%.3g %.3g %.3g" % (np.var(u), np.var(v), np.var(rho)))
        path = "cyl.%09d.raw" % cycle
        vort = np.roll(u, [0, 1], [0, 1]) - np.roll(
            u, [0, -1], [0, 1]) - np.roll(v, [1, 0], [0, 1]) + np.roll(
                v, [-1, 0], [0, 1])
        with open(path, "wb") as fid:
            for field in [u, v, rho, vort]:
                field[xb, yb] = np.nan
                fid.write(field.tobytes("F"))
