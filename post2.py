import matplotlib.pyplot as plt
import matplotlib
import sys
import numpy as np
import re

nx, ny = 64 * 32, 128 * 32
dtype = np.float64

for path in sys.argv[1:]:
    u, v, rho, vort = np.memmap(path, dtype=dtype, shape=(4, nx, ny))
    plt.axis("off")
    plt.axis("equal")
    plt.tight_layout()
    lo, hi = np.quantile(vort, [0.05, 0.95])
    plt.imshow(vort, matplotlib.cm.jet, vmin=lo, vmax=hi)
    png = re.sub("\.raw$", "", path) + ".png"
    print(png)
    plt.savefig(png, bbox_inches=0, dpi=500)
    plt.close()
