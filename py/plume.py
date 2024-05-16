import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import ticker

import seaborn as sb
import sys

def read_dump(fname):
    count = 0; lim = 9;
    arr = [];
    f = open(fname);
    for line in f: 
        if(count >= lim):
            arr.append([float(x) for x in line.split()])
        count =  count + 1;
    arr = np.asarray(arr);
    f.close()
    return arr

fname = sys.argv[1]
outDSMC = read_dump(fname);
outDSMC = outDSMC[outDSMC[:, 0].argsort()]

xDSMC   = outDSMC[:, 0]
yDSMC   = outDSMC[:, 1]
zDSMC   = outDSMC[:, 2]

NDSMC   = outDSMC[:, 3]
nDSMC   = outDSMC[:, 4]
uDSMC   = outDSMC[:, 5]
vDSMC   = outDSMC[:, 6]
wDSMC   = outDSMC[:, 7]
TDSMC   = outDSMC[:, 8]
pDSMC   = outDSMC[:, 9]
txyDSMC = outDSMC[:,10]
tyzDSMC = outDSMC[:,11]
txzDSMC = outDSMC[:,12]
qxDSMC  = outDSMC[:,13]
qyDSMC  = outDSMC[:,14]
qzDSMC  = outDSMC[:,15]

nx = 100
ny = 100
nz = 40

rnz = 1e-3;
L = 40*rnz;
Lx = 2.5*L;
Ly = 2.5*L;
Lz = L;

dx = Lx/nx;
dy = Ly/ny;
dz = Lz/nz;


yL = 0;
yH = dy;
indy = np.array(np.where(np.logical_and(yDSMC>yL, yDSMC<yH)));
indy = indy.flatten()

nslice = np.zeros((nx,nz))
wslice = np.zeros((nx,nz))
for ix in range(nx):
    xL = ix*dx;
    xH = (ix+1)*dx;
    indx = np.array(np.where(np.logical_and(xDSMC>xL, xDSMC<xH)));
    indx = indx.flatten()

    for iz in range(nz):
        zL = iz*dz;
        zH = (iz+1)*dz;
        indz = np.array(np.where(np.logical_and(zDSMC>zL, zDSMC<zH)));
        indz = indz.flatten()
    
        ind = np.intersect1d(indx,indz)
        ind = np.intersect1d(ind,indy)
        nslice[ix,iz] = np.mean(NDSMC[ind])
        wslice[ix,iz] = np.mean(wDSMC[ind])

#############################################################
fig,(ax0,ax1) = plt.subplots(nrows=1,ncols=2);

#############################################################
sb.heatmap(-wslice,cmap="jet",ax=ax0,vmin=-130,vmax=740)
nmax = np.max(np.log(nslice))
nmin = np.min(np.log(nslice))
sb.heatmap(np.log(nslice),cmap="jet",ax=ax1,vmin=nmin,vmax=nmax)

fig.set_size_inches(8.5, 8.5)

mng = plt.get_current_fig_manager()
mng.full_screen_toggle()
plt.savefig(sys.argv[2],bbox_inches='tight')

sys.exit()









