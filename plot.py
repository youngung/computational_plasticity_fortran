import numpy as np
import matplotlib.pyplot as plt
fig=plt.figure()
ax=fig.add_subplot(111)

dat=np.loadtxt('ys.txt').T
ax.plot(dat[0],dat[1],label=r'$t\sigma$')
ax.plot(dat[3],dat[4],label=r'$r\sigma$ (Hill)')
ax.plot(dat[6],dat[7],'--',label=r'$r\sigma$ (von Mises)')

ax.legend()
ax.set_aspect('equal')
ax.grid()
fig.savefig('ys.png')
