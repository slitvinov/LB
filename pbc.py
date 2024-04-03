import torch as np
if np.cuda.is_available():
    np.set_default_device("cuda")

nx = 800
ny = 200
nsteps = 10000
nplot = 1000
x0 = nx // 5
y0 = ny // 2
r0 = ny // 20
u0 = 0.1
Re = 100
nu = u0 * 2. * r0 / Re
omega = 1. / (3 * nu + 1. / 2.)
t = 4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36
cx = 0, 1, 0, -1, 0, 1, -1, -1, 1
cy = 0, 0, 1, 0, -1, 1, 1, -1, -1
x, y = np.meshgrid(np.arange(nx + 2), np.arange(ny + 2), indexing="ij")
obst = (x - x0)**2 + (y - y0)**2 <= r0**2
u = np.zeros((nx + 2, ny + 2), dtype=np.float64)
v = np.zeros((nx + 2, ny + 2), dtype=np.float64)
rho0 = 1
rho = np.full((nx + 2, ny + 2), rho0, dtype=np.float64)
f = [np.empty((nx + 2, ny + 2), dtype=np.float64) for i in range(9)]
feq = np.empty((nx + 2, ny + 2), dtype=np.float64)
c = np.empty((nx + 2, ny + 2), dtype=np.float64)

for i in range(9):
    c[:] = 3 * (cx[i] * u + cy[i] * v)
    f[i][:] = rho0 * t[i] * (1 + c + 1 / 2 * (c * c) - 3 / 2 * (u**2 + v**2))

for istep in range(nsteps):
    f[1][0, 1:ny + 1] = f[1][nx, 1:ny + 1]
    f[5][0, 1:ny + 1] = f[5][nx, 1:ny + 1]
    f[8][0, 1:ny + 1] = f[8][nx, 1:ny + 1]

    f[3][nx + 1, 1:ny + 1] = f[3][1, 1:ny + 1]
    f[6][nx + 1, 1:ny + 1] = f[6][1, 1:ny + 1]
    f[7][nx + 1, 1:ny + 1] = f[7][1, 1:ny + 1]

    f[2][1:nx + 1, 0] = f[2][1:nx + 1, ny]
    f[5][1:nx + 1, 0] = f[5][1:nx + 1, ny]
    f[6][1:nx + 1, 0] = f[6][1:nx + 1, ny]

    f[4][1:nx + 1, ny + 1] = f[4][1:nx + 1, 1]
    f[7][1:nx + 1, ny + 1] = f[7][1:nx + 1, 1]
    f[8][1:nx + 1, ny + 1] = f[8][1:nx + 1, 1]


    rho[:] = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8]
    u[:] = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) / rho
    v[:] = (f[5] + f[2] + f[6] - f[7] - f[4] - f[8]) / rho

    force = 1e-5
    f[1][1:nx+1,1:ny+1] += force
    f[5][1:nx+1,1:ny+1] += force
    f[8][1:nx+1,1:ny+1] += force

    f[3][1:nx+1,1:ny+1] -= force
    f[6][1:nx+1,1:ny+1] -= force
    f[7][1:nx+1,1:ny+1] -= force

    # collision and bounce-back boundary condition
    for i in range(9):
        c[:] = 3 * (cx[i] * u + cy[i] * v)
        feq[:] = rho * t[i] * (1 + c + 1 / 2 * (c * c) - 3 / 2 * (u**2 + v**2))
        f[i][~obst] = f[i][~obst] * (1 - omega) + omega * feq[~obst]
    f[1][obst], f[3][obst] = f[3][obst], f[1][obst]
    f[2][obst], f[4][obst] = f[4][obst], f[2][obst]
    f[5][obst], f[7][obst] = f[7][obst], f[5][obst]
    f[6][obst], f[8][obst] = f[8][obst], f[6][obst]

    # streaming step
    for i in range(9):
        f[i][:] = np.roll(f[i], (cx[i], cy[i]), (0, 1))

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
                field[obst] = np.nan
                if hasattr(field, "numpy"):
                    field = field.cpu().detach().numpy()
                fid.write(field[1:nx+1, 1:ny+1].tobytes("F"))
"""

625
301
748

"""

