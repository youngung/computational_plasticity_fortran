import numpy as np
import matplotlib.pyplot as plt
fig=plt.figure(figsize=(3,3))
ax=fig.add_subplot(111)

dat=np.loadtxt('et.txt').T
ax.plot(dat[0],dat[1],'-',label=r'$\varepsilon(t)$ (analytical)',zorder=10)

dat=np.loadtxt('et_num.txt').T
ax.plot(dat[0],dat[1],'o',mfc='None',label=r'$\varepsilon(t)$ (numerical)')

ax.set_xlabel('time [sec]')
ax.set_ylabel(r'$\varepsilon$')
ax.legend()
#ax.set_aspect('equal')
plt.tight_layout()
ax.grid()
fig.savefig('x.png')
