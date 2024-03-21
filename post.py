import math
import matplotlib.animation
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def update(path):
    size = os.path.getsize(path)
    ny = math.isqrt(size // dtype.itemsize // nfields // 2)
    nx = 2 * ny
    u, v, rho, vort = np.memmap(path, dtype=dtype, shape=(4, ny, nx))
    plt.imshow(u, matplotlib.cm.jet)
    plt.tight_layout()

dtype = np.dtype("float64")
nfields = 4
plt.axis("off")
plt.axis("equal")
ani = matplotlib.animation.FuncAnimation(fig=plt.gcf(),
                                         func=update,
                                         frames=sys.argv[1:])
ani.save("u.mp4")
