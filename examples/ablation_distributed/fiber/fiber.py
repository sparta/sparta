import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import ticker

import sys

plt.rc('text', usetex=False)
plt.rc('font', size=25)
plt.rc('axes', titlesize=20)
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['lines.markersize'] = 6.0
plt.rcParams['lines.markeredgecolor'] = 'k'
plt.rcParams['lines.markeredgewidth'] = 1.0
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['legend.title_fontsize'] = 'x-small'

def datread(fname):
    count = 0;
    arr = [];
    f = open(fname, 'r');
    for line in f:
      if(count > 0):
        arr.append([float(x) for x in line.split()])
      count = count + 1
    arr = np.asarray(arr);
    f.close()
    return arr

prefix_rxn = sys.argv[1]
prefix_diff = sys.argv[2]
outname = sys.argv[3]

# volume

aoR = datread(prefix_rxn+'case0/vol.dat')
ioR = datread(prefix_rxn+'case1/vol.dat')
amR = datread(prefix_rxn+'case2/vol.dat')
imR = datread(prefix_rxn+'case3/vol.dat')

aoD = datread(prefix_diff+'case0/vol.dat')
ioD = datread(prefix_diff+'case1/vol.dat')
amD = datread(prefix_diff+'case2/vol.dat')
imD = datread(prefix_diff+'case3/vol.dat')

A = (12.5e-6+2.5e-6)*(12.5e-6+2.5e-6)*42.0e-6

fig,ax = plt.subplots(1,2)

# rxn

sz = len(aoR)
x = np.arange(sz)
ax[0].plot(x,(A-aoR)/(A-aoR[0]),'-k',linewidth=3,label=r'$Ave/Single$')

sz = len(ioR)
x = np.arange(sz)
ax[0].plot(x,(A-ioR)/(A-ioR[0]),'-r',linewidth=3,label=r'$Inner/Single$')

sz = len(amR)
x = np.arange(sz)
ax[0].plot(x,(A-amR)/(A-amR[0]),'-g',linewidth=3,label=r'$Ave/Multi$')

sz = len(imR)
x = np.arange(sz)
ax[0].plot(x,(A-imR)/(A-imR[0]),'-b',linewidth=3,label=r'$Inner/Multi$')

# diff

sz = len(aoD)
x = np.arange(sz)
ax[1].plot(x,(A-aoD)/(A-aoD[0]),'-k',linewidth=3,label=r'$Ave/Single$')

sz = len(ioD)
x = np.arange(sz)
ax[1].plot(x,(A-ioD)/(A-ioD[0]),'-r',linewidth=3,label=r'$Inner/Single$')

sz = len(amD)
x = np.arange(sz)
ax[1].plot(x,(A-amD)/(A-amD[0]),'-g',linewidth=3,label=r'$Ave/Multi$')

sz = len(imD)
x = np.arange(sz)
ax[1].plot(x,(A-imD)/(A-imD[0]),'-b',linewidth=3,label=r'$Inner/Multi$')

# format

#ax[0].set_xlabel(r'$N_{ablate}$')
#ax[1].set_xlabel(r'$N_{ablate}$')
ax[0].set_ylabel(r'$m/m_0$')

ax[0].set_ylim([0.0,1.0])
ax[0].set_ylim([0.4,1.0])

for i in range(2):
  ax[i].set_xlim([0,600])

  ax[i].minorticks_on()
  ax[i].tick_params(axis='both',direction='in',which='both',right=True,top=True)

  xval, yval = ax[i].get_xlim(), ax[i].get_ylim()
  xlims = xval[1]-xval[0];
  ylims = yval[1]-yval[0];
  ax[i].set_aspect(1.0*(xlims/ylims), adjustable='box')

ax[0].legend(prop={'size': 16})
fig.set_size_inches(12.,12.)

fig.tight_layout();
mng = plt.get_current_fig_manager()
mng.full_screen_toggle()
plt.savefig(outname)



