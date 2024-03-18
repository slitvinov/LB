import matplotlib.pyplot as plt
import matplotlib
import matplotlib.animation
import sys
import numpy as np


def update(path):
    dtype = np.float64
    u, v, rho, vort = np.memmap(path, dtype=dtype, shape=(4, nx, ny))
    plt.imshow(u, matplotlib.cm.jet)
    plt.tight_layout()


nx, ny = 64, 128
plt.axis("off")
plt.axis("equal")
ani = matplotlib.animation.FuncAnimation(fig=plt.gcf(),
                                         func=update,
                                         frames=sys.argv[1:])
ani.save("u.mp4")
