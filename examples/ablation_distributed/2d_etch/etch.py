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

def fread(fname):
    count = 0;
    pts = [];
    lines = [];
    npts = 0
    nlines = 0
    f = open(fname, 'r');
    for line in f:
      if(count == 2):
        npts = int(line.split()[0])
      elif(count == 3):
        nlines = int(line.split()[0])
      elif(count >= 7 and count <= (6+npts)):
        pts.append([float(x) for x in line.split()])
      elif(count > (9+npts)):
        lines.append([int(x) for x in line.split()])
      count = count + 1
    pts = np.asarray(pts);
    lines = np.asarray(lines);
    f.close()
    return pts,lines

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

prefix = sys.argv[1] # prefix
outname0 = sys.argv[2]

# surfaces

ao0pts,ao0l = fread(prefix+'etch.0.ave.one.surf')
ao1pts,ao1l = fread(prefix+'etch.1.ave.one.surf')
ao2pts,ao2l = fread(prefix+'etch.2.ave.one.surf')
ao3pts,ao3l = fread(prefix+'etch.3.ave.one.surf')
ao4pts,ao4l = fread(prefix+'etch.4.ave.one.surf')
ao5pts,ao5l = fread(prefix+'etch.5.ave.one.surf')

io0pts,io0l = fread(prefix+'etch.0.inner.one.surf')
io1pts,io1l = fread(prefix+'etch.1.inner.one.surf')
io2pts,io2l = fread(prefix+'etch.2.inner.one.surf')
io3pts,io3l = fread(prefix+'etch.3.inner.one.surf')
io4pts,io4l = fread(prefix+'etch.4.inner.one.surf')
io5pts,io5l = fread(prefix+'etch.5.inner.one.surf')

am0pts,am0l = fread(prefix+'etch.0.ave.multi.surf')
am1pts,am1l = fread(prefix+'etch.1.ave.multi.surf')
am2pts,am2l = fread(prefix+'etch.2.ave.multi.surf')
am3pts,am3l = fread(prefix+'etch.3.ave.multi.surf')
am4pts,am4l = fread(prefix+'etch.4.ave.multi.surf')
am5pts,am5l = fread(prefix+'etch.5.ave.multi.surf')

im0pts,im0l = fread(prefix+'etch.0.inner.multi.surf')
im1pts,im1l = fread(prefix+'etch.1.inner.multi.surf')
im2pts,im2l = fread(prefix+'etch.2.inner.multi.surf')
im3pts,im3l = fread(prefix+'etch.3.inner.multi.surf')
im4pts,im4l = fread(prefix+'etch.4.inner.multi.surf')
im5pts,im5l = fread(prefix+'etch.5.inner.multi.surf')

# plot

# colors
r = np.linspace(0,1,6)
g = np.linspace(0,1,6)
b = np.linspace(0,1,6)

g = g*0
r = 1.0 - b;

c0 = (r[0],g[0],b[0])
c1 = (r[1],g[1],b[1])
c2 = (r[2],g[2],b[2])
c3 = (r[3],g[3],b[3])
c4 = (r[4],g[4],b[4])
c5 = (r[5],g[5],b[5])

fig,ax = plt.subplots(2,2)

############################################################################
# ave - single

sz0 = np.shape(ao0l)[0]
sz1 = np.shape(ao1l)[0]
sz2 = np.shape(ao2l)[0]
sz3 = np.shape(ao3l)[0]
sz4 = np.shape(ao4l)[0]
sz5 = np.shape(ao5l)[0]
szmax = np.max([sz0,sz1,sz2,sz3,sz4,sz5])

for i in range(szmax):
  if(i < sz0):
    p1 = ao0l[i,1]-1  
    p2 = ao0l[i,2]-1
    if(i==0):
      ax[0][0].plot([ao0pts[p1,1],ao0pts[p2,1]],
                    [ao0pts[p1,2],ao0pts[p2,2]],
                    '-',color=c0,label=r'$0$')
    else:
      ax[0][0].plot([ao0pts[p1,1],ao0pts[p2,1]],
                    [ao0pts[p1,2],ao0pts[p2,2]],
                    '-',color=c0)

  if(i < sz1):
    p1 = ao1l[i,1]-1
    p2 = ao1l[i,2]-1
    if(i==0):
      ax[0][0].plot([ao1pts[p1,1],ao1pts[p2,1]],
                    [ao1pts[p1,2],ao1pts[p2,2]],
                     '-',color=c1,label=r'$125$')
    else:
      ax[0][0].plot([ao1pts[p1,1],ao1pts[p2,1]],
                    [ao1pts[p1,2],ao1pts[p2,2]],
                     '-',color=c1)

  if(i < sz2):
    p1 = ao2l[i,1]-1
    p2 = ao2l[i,2]-1
    if(i==0):
      ax[0][0].plot([ao2pts[p1,1],ao2pts[p2,1]],
                    [ao2pts[p1,2],ao2pts[p2,2]],
                     '-',color=c2,label=r'$250$')
    else:
      ax[0][0].plot([ao2pts[p1,1],ao2pts[p2,1]],
                    [ao2pts[p1,2],ao2pts[p2,2]],
                     '-',color=c2)

  if(i < sz3):
    p1 = ao3l[i,1]-1
    p2 = ao3l[i,2]-1
    if(i==0):
      ax[0][0].plot([ao3pts[p1,1],ao3pts[p2,1]],
                    [ao3pts[p1,2],ao3pts[p2,2]],
                     '-',color=c3,label=r'$375$')
    else:
      ax[0][0].plot([ao3pts[p1,1],ao3pts[p2,1]],
                    [ao3pts[p1,2],ao3pts[p2,2]],
                     '-',color=c3)

  if(i < sz4):
    p1 = ao4l[i,1]-1
    p2 = ao4l[i,2]-1
    if(i==0):
      ax[0][0].plot([ao4pts[p1,1],ao4pts[p2,1]],
                    [ao4pts[p1,2],ao4pts[p2,2]],
                     '-',color=c4,label=r'$500$')
    else:
      ax[0][0].plot([ao4pts[p1,1],ao4pts[p2,1]],
                    [ao4pts[p1,2],ao4pts[p2,2]],
                     '-',color=c4)

  if(i < sz5):
    p1 = ao5l[i,1]-1
    p2 = ao5l[i,2]-1
    if(i==0):
      ax[0][0].plot([ao5pts[p1,1],ao5pts[p2,1]],
                    [ao5pts[p1,2],ao5pts[p2,2]],
                     '-',color=c5,label=r'$625$')
    else:
      ax[0][0].plot([ao5pts[p1,1],ao5pts[p2,1]],
                    [ao5pts[p1,2],ao5pts[p2,2]],
                     '-',color=c5)

############################################################################
# inner - single

sz0 = np.shape(io0l)[0]
sz1 = np.shape(io1l)[0]
sz2 = np.shape(io2l)[0]
sz3 = np.shape(io3l)[0]
sz4 = np.shape(io4l)[0]
sz5 = np.shape(io5l)[0]
szmax = np.max([sz0,sz1,sz2,sz3,sz4,sz5])

for i in range(szmax):
  if(i < sz0):
    p1 = io0l[i,1]-1  
    p2 = io0l[i,2]-1
    ax[0][1].plot([io0pts[p1,1],io0pts[p2,1]],
                  [io0pts[p1,2],io0pts[p2,2]],
                   '-',color=c0)

  if(i < sz1):
    p1 = io1l[i,1]-1
    p2 = io1l[i,2]-1
    ax[0][1].plot([io1pts[p1,1],io1pts[p2,1]],
                  [io1pts[p1,2],io1pts[p2,2]],
                   '-',color=c1)

  if(i < sz2):
    p1 = io2l[i,1]-1
    p2 = io2l[i,2]-1
    ax[0][1].plot([io2pts[p1,1],io2pts[p2,1]],
                  [io2pts[p1,2],io2pts[p2,2]],
                   '-',color=c2)

  if(i < sz3):
    p1 = io3l[i,1]-1
    p2 = io3l[i,2]-1
    ax[0][1].plot([io3pts[p1,1],io3pts[p2,1]],
                  [io3pts[p1,2],io3pts[p2,2]],
                   '-',color=c3)

  if(i < sz4):
    p1 = io4l[i,1]-1
    p2 = io4l[i,2]-1
    ax[0][1].plot([io4pts[p1,1],io4pts[p2,1]],
                  [io4pts[p1,2],io4pts[p2,2]],
                   '-',color=c4)

  if(i < sz5):
    p1 = io5l[i,1]-1
    p2 = io5l[i,2]-1
    ax[0][1].plot([io5pts[p1,1],io5pts[p2,1]],
                  [io5pts[p1,2],io5pts[p2,2]],
                   '-',color=c5)

############################################################################
# ave - multiple

sz0 = np.shape(am0l)[0]
sz1 = np.shape(am1l)[0]
sz2 = np.shape(am2l)[0]
sz3 = np.shape(am3l)[0]
sz4 = np.shape(am4l)[0]
sz5 = np.shape(am5l)[0]
szmax = np.max([sz0,sz1,sz2,sz3,sz4,sz5])

for i in range(szmax):
  if(i < sz0):
    p1 = am0l[i,1]-1  
    p2 = am0l[i,2]-1
    ax[1][0].plot([am0pts[p1,1],am0pts[p2,1]],
                  [am0pts[p1,2],am0pts[p2,2]],
                   '-',color=c0)

  if(i < sz1):
    p1 = am1l[i,1]-1
    p2 = am1l[i,2]-1
    ax[1][0].plot([am1pts[p1,1],am1pts[p2,1]],
                  [am1pts[p1,2],am1pts[p2,2]],
                   '-',color=c1)

  if(i < sz2):
    p1 = am2l[i,1]-1
    p2 = am2l[i,2]-1
    ax[1][0].plot([am2pts[p1,1],am2pts[p2,1]],
                  [am2pts[p1,2],am2pts[p2,2]],
                   '-',color=c2)

  if(i < sz3):
    p1 = am3l[i,1]-1
    p2 = am3l[i,2]-1
    ax[1][0].plot([am3pts[p1,1],am3pts[p2,1]],
                  [am3pts[p1,2],am3pts[p2,2]],
                   '-',color=c3)

  if(i < sz4):
    p1 = am4l[i,1]-1
    p2 = am4l[i,2]-1
    ax[1][0].plot([am4pts[p1,1],am4pts[p2,1]],
                  [am4pts[p1,2],am4pts[p2,2]],
                   '-',color=c4)

  if(i < sz5):
    p1 = am5l[i,1]-1
    p2 = am5l[i,2]-1
    ax[1][0].plot([am5pts[p1,1],am5pts[p2,1]],
                  [am5pts[p1,2],am5pts[p2,2]],
                   '-',color=c5)

############################################################################
# inner - multiple

sz0 = np.shape(im0l)[0]
sz1 = np.shape(im1l)[0]
sz2 = np.shape(im2l)[0]
sz3 = np.shape(im3l)[0]
sz4 = np.shape(im4l)[0]
sz5 = np.shape(im5l)[0]
szmax = np.max([sz0,sz1,sz2,sz3,sz4,sz5])

for i in range(szmax):
  if(i < sz0):
    p1 = im0l[i,1]-1  
    p2 = im0l[i,2]-1
    ax[1][1].plot([im0pts[p1,1],im0pts[p2,1]],
                  [im0pts[p1,2],im0pts[p2,2]],
                   '-',color=c0)

  if(i < sz1):
    p1 = im1l[i,1]-1
    p2 = im1l[i,2]-1
    ax[1][1].plot([im1pts[p1,1],im1pts[p2,1]],
                  [im1pts[p1,2],im1pts[p2,2]],
                   '-',color=c1)

  if(i < sz2):
    p1 = im2l[i,1]-1
    p2 = im2l[i,2]-1
    ax[1][1].plot([im2pts[p1,1],im2pts[p2,1]],
                  [im2pts[p1,2],im2pts[p2,2]],
                   '-',color=c2)

  if(i < sz3):
    p1 = im3l[i,1]-1
    p2 = im3l[i,2]-1
    ax[1][1].plot([im3pts[p1,1],im3pts[p2,1]],
                  [im3pts[p1,2],im3pts[p2,2]],
                   '-',color=c3)

  if(i < sz4):
    p1 = im4l[i,1]-1
    p2 = im4l[i,2]-1
    ax[1][1].plot([im4pts[p1,1],im4pts[p2,1]],
                  [im4pts[p1,2],im4pts[p2,2]],
                   '-',color=c4)

  if(i < sz5):
    p1 = im5l[i,1]-1
    p2 = im5l[i,2]-1
    ax[1][1].plot([im5pts[p1,1],im5pts[p2,1]],
                  [im5pts[p1,2],im5pts[p2,2]],
                   '-',color=c5)

############################################################################

ax[0][0].set_title('AS')
ax[0][1].set_title('IS')
ax[1][0].set_title('AM')
ax[1][1].set_title('IM')

for i in range(2):
  for j in range(2):
    ax[i][j].plot([11,11],[0,10],'-k',linewidth=1)
    ax[i][j].plot([0,22],[8,8],':k',linewidth=1)
    ax[i][j].plot([0,22],[7,7],':k',linewidth=1)
    ax[i][j].plot([0,22],[6,6],':k',linewidth=1)
    ax[i][j].plot([0,22],[5,5],':k',linewidth=1)
    ax[i][j].plot([0,22],[4,4],':k',linewidth=1)
    #ax[i][j].plot([0,22],[3,3],':k',linewidth=1)
    ax[i][j].set_ylim([2.5,8.5])
    ax[i][j].set_xlim([0,21.5])
    ax[i][j].set_xticks([])
    ax[i][j].set_yticks([3,4,5,6,7,8])

ax[0][0].set_yticklabels([3,4,5,6,7,8])
ax[1][0].set_yticklabels([3,4,5,6,7,8])

ax[0][0].set_ylabel(r'$y$')
ax[1][0].set_ylabel(r'$y$')
ax[1][0].set_xlabel(r'$x$')
ax[1][1].set_xlabel(r'$x$')

'''
ax.set_xlabel(r'$t\;[ms]$')
ax.set_ylabel(r'$m/m_0$')

ax.set_xlim([0,500])
ax.set_ylim([1.0,1.25])
ax.set_xticks([0,100,200,300,400,500])
ax.set_yticks([1.0,1.05,1.10,1.15,1.20,1.25])
ax.tick_params(right=True, top=True)

xval, yval = ax.get_xlim(), ax.get_ylim()
xlims = xval[1]-xval[0];
ylims = yval[1]-yval[0];
ax.set_aspect(1.0*(xlims/ylims), adjustable='box')
'''
ax[0][0].legend(title=r'$N_{ablate}$',ncol=2)
fig.set_size_inches(12.,12.)

fig.tight_layout();
mng = plt.get_current_fig_manager()
mng.full_screen_toggle()
plt.savefig(outname0)
