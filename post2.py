import matplotlib.pyplot as plt
import matplotlib
import sys
import numpy as np
import re
import os
import math

dtype = np.dtype("float64")
nfields = 4
for path in sys.argv[1:]:
    size = os.path.getsize(path)
    ny = math.isqrt(size // dtype.itemsize // nfields // 4)
    nx = 4 * ny
    u, v, rho, vort = np.memmap(path, dtype=dtype, shape=(4, ny, nx))
    plt.axis("off")
    plt.axis("equal")
    plt.tight_layout()
    lo, hi = np.quantile(u, [0.05, 0.95])
    plt.imshow(u, matplotlib.cm.jet, vmin=lo, vmax=hi)
    png = re.sub("\.raw$", "", path) + ".png"
    print(png)
    plt.savefig(png, bbox_inches=0, dpi=500)
    plt.close()
