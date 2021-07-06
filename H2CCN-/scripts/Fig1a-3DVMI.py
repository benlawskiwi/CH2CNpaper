from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

print("Fig1: warning takes a little time to run ...")
IM = np.loadtxt("../data/CH2CN-6N3CMb512.txt")

n, m = IM.shape
n /= 2
m /= 2

Xr = np.arange(-m, m, 1)
Yr = np.arange(-n, n, 1)
X, Y = np.meshgrid(Xr, Yr)

IM = IM.T   # to give the correct orientation

width = 3.33
height = width   #  equal
one_col = (width, height)

fig = plt.figure(figsize=one_col)
ax = fig.gca(projection='3d')

ax.view_init(elev=70, azim=45) # 70, 40)

surf = ax.plot_surface(X, Y, IM, rstride=1, cstride=1, cmap=cm.jet, #cm.viridis,
                       linewidth=0, antialiased=False, shade=False,
                       facecolor='w', vmin=0)
#ax.text(200, 0, 0, r'$1^0_0 2^0_0$', (0,1,1) , color='y', fontsize='larger')

#ax._axis3don = True
ax.grid(False)
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

plt.tick_params(axis='both', left='off', right='off')

ax.axis(xmin=-m, xmax=m, ymin=-n, ymax=n, zmin=0)
ax.set_zticks([])
ax.w_zaxis.line.set_lw(0)

# colour bar
cf = fig.add_axes([0.17, 0.47, 0.03, 0.25])
cbar = fig.colorbar(surf, ticks=[0, 250], cax=cf)
cf.yaxis.set_ticks_position('right')
cf.yaxis.set_label_position('left')
cf.set_ylabel('electrons')


ax.set_xticks([-n, 0, n])
ax.set_xticklabels([])
ax.set_xlabel(r"laser polarization $\longrightarrow$", labelpad=3, ha='right')

ax.set_ylabel(r"radius (pixels)", labelpad=5)
ax.set_yticks([-m, 0, m])

ax.set_zlim([0, 250])

fig.set_size_inches(*one_col)

plt.savefig("Fig1a-3DVMI.tif", dpi=300, bb_inches='tight')

plt.show()
