import matplotlib.pyplot as plt
import matplotlib
import matplotlib.animation
import sys
import numpy as np

def update(path):
    u = np.memmap(path, dtype=np.float64, shape=(nx, ny), order='C')
    plt.imshow(u,  matplotlib.cm.jet)

nx, ny = 64, 128
plt.axis("off")
plt.tight_layout()
ani = matplotlib.animation.FuncAnimation(fig=plt.gcf(), func=update, frames=sys.argv[1:])
ani.save("u.gif")
