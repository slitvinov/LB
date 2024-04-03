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
    n = math.isqrt(size // dtype.itemsize // nfields)
    u, v, rho, vort = np.memmap(path, dtype=dtype, shape=(nfields, n, n))
    plt.figure(figsize=(20, 5), dpi=150, frameon=False)
    plt.axis("off")
    plt.axis("equal")
    plt.tight_layout()
    plt.imshow(u, matplotlib.cm.jet)
    png = re.sub("\.raw$", "", path) + ".png"
    print(png)
    plt.savefig(png, bbox_inches=0)
    plt.close()
