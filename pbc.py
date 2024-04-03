import math
import torch as np
if np.cuda.is_available():
    np.set_default_device("cuda")

n = 1000
nsteps = 10000
nplot = 100
cs2 = 1 / 3
omega = 3 / 2
t = [4 / 9] + [1 / 9] * 4 + [1 / 36] * 4
cx = 0, 1, 0, -1, 0, 1, -1, -1, 1
cy = 0, 0, 1, 0, -1, 1, 1, -1, -1
x = np.asarray([ 2 * math.pi * (i - 1) / n for i in range(n + 2) ], dtype=np.float64)
u = np.empty((n + 2, n + 2), dtype=np.float64)
v = np.empty((n + 2, n + 2), dtype=np.float64)
u0 = 0.1
visc = (1 / omega - 0.5) * cs2
print(visc, re)
np.outer(u0 * np.sin(x), np.cos(x), out=u)
np.outer(-u0 * np.cos(x), np.sin(x), out=v)
rho = np.ones((n + 2, n + 2), dtype=np.float64)
f = [np.empty((n + 2, n + 2), dtype=np.float64) for i in range(9)]
feq = np.empty((n + 2, n + 2), dtype=np.float64)
c = np.empty((n + 2, n + 2), dtype=np.float64)
for i in range(9):
    c[:] = (cx[i] * u + cy[i] * v) / cs2
    f[i][:] = t[i] * (1 + c + 1 / 2 * (c * c) - 3 / 2 * (u**2 + v**2))

for istep in range(nsteps):
    f[1][0, 1:n + 1] = f[1][n, 1:n + 1]
    f[5][0, 1:n + 1] = f[5][n, 1:n + 1]
    f[8][0, 1:n + 1] = f[8][n, 1:n + 1]

    f[3][n + 1, 1:n + 1] = f[3][1, 1:n + 1]
    f[6][n + 1, 1:n + 1] = f[6][1, 1:n + 1]
    f[7][n + 1, 1:n + 1] = f[7][1, 1:n + 1]

    f[2][1:n + 1, 0] = f[2][1:n + 1, n]
    f[5][1:n + 1, 0] = f[5][1:n + 1, n]
    f[6][1:n + 1, 0] = f[6][1:n + 1, n]

    f[4][1:n + 1, n + 1] = f[4][1:n + 1, 1]
    f[7][1:n + 1, n + 1] = f[7][1:n + 1, 1]
    f[8][1:n + 1, n + 1] = f[8][1:n + 1, 1]

    f[5][0, 0] = f[5][n, n]
    f[7][n + 1, n + 1] = f[7][1, 1]

    f[5][0, 0] = f[5][n, n]
    f[7][n + 1, n + 1] = f[7][1, 1]

    f[8][n + 1, 0] = f[8][1, n]
    f[6][0, n + 1] = f[6][n, 1]

    rho[:] = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8]
    u[:] = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) / rho
    v[:] = (f[5] + f[2] + f[6] - f[7] - f[4] - f[8]) / rho

    # collision and bounce-back boundary condition
    for i in range(9):
        c[:] = (cx[i] * u + cy[i] * v) / cs2
        feq[:] = rho * t[i] * (1 + c + 1 / 2 * (c * c) - 3 / 2 * (u**2 + v**2))
        f[i][:] = f[i] * (1 - omega) + omega * feq

    # streaming step
    for i in range(9):
        f[i][:][1:n + 1, 1:n + 1] = np.roll(f[i], (cx[i], cy[i]),
                                              (0, 1))[1:n + 1, 1:n + 1]

    if istep % nplot == 0:
        path = "cyl.%09d.raw" % istep
        print("cylinder.py: %s" % path)
        vort = np.roll(u, [0, 1], [0, 1]) - np.roll(
            u, [0, -1], [0, 1]) - np.roll(v, [1, 0], [0, 1]) + np.roll(
                v, [-1, 0], [0, 1])
        for name, field in ["u", u], ["v", v], ["rho", rho], ["vort", vort]:
            print(
                "cylinder.py: %10s: mean,min,max,std: %+.3e %+.3e %+.3e %+.3e"
                % (name, np.mean(field), np.min(field), np.max(field),
                   np.std(field)))
        with open(path, "wb") as fid:
            for field in [u, v, rho, vort]:
                if hasattr(field, "numpy"):
                    field = field.cpu().detach().numpy()
                fid.write(field[1:n + 1, 1:n + 1].tobytes("F"))
"""

625
301
748

oooo
oxxo
oxxo
oooo

"""
