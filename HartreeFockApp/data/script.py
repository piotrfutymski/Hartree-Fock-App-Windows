import matplotlib.pyplot as plt
import numpy as np
import sys

# generate 2 2d grids for the x & y bounds

for a in range(1, len(sys.argv)):
    y, x = np.meshgrid(np.linspace(-3.0, 3.0, 601), np.linspace(-3.0, 3.0, 601))
    z = x*y
    f = open(sys.argv[a]+".txt", "r")
    for i in range(0, 601):
        for j in range(0,601):
            z[i,j] = f.readline()
    f.close()
    z = z[:-1, :-1]
    z_min, z_max = -np.abs(z).max(), np.abs(z).max()

    fig, ax = plt.subplots()

    c = ax.pcolormesh(x, y, z, cmap='RdBu', vmin=z_min, vmax=z_max)
    ax.set_title(sys.argv[a])

    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(c, ax=ax)

    plt.savefig("pic/"+sys.argv[a]+".png")



