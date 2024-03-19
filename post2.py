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
    ny = math.isqrt(size // dtype.itemsize // nfields // 2)
    nx = 2 * ny
    u, v, rho, vort = np.memmap(path, dtype=dtype, shape=(4, ny, nx))
    plt.axis("off")
    plt.axis("equal")
    plt.tight_layout()
    plt.imshow(u, matplotlib.cm.jet)
    png = re.sub("\.raw$", "", path) + ".png"
    print(png)
    plt.savefig(png, bbox_inches=0, dpi=500)
    plt.close()
