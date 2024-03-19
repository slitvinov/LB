import matplotlib.pyplot as plt
import matplotlib
import sys
import numpy as np
import re

nx, ny = 64 * 8, 128 * 8
plt.axis("off")
plt.axis("equal")
plt.tight_layout()
dtype = np.float64

for path in sys.argv[1:]:
    u, v, rho, vort = np.memmap(path, dtype=dtype, shape=(4, nx, ny))
    plt.imshow(u, matplotlib.cm.jet)
    png = re.sub("\.raw$", "", path) + ".png"
    plt.savefig(png)
