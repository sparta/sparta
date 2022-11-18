#!/usr/bin/env python

# Script:  grid_refine.py
# Purpose: create a refined hierarchical grid around
#          a SPARTA surf file, based on closeness to surf
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# Syntax:  grid_refine.py switch args ...
#          -s surffile = SPARTA surf file (no default)
#          -b xlo xhi ylo yhi zlo zhi = simulation box bounds
#             (def = 0 1 0 1 0 1)
#             zlo/zhi must be specified, but are ignored for 2d
#          -n Nx Ny Nz = initial coarse grid (def = 2 2 2)
#                        Nz set to 1 for 2d
#          -m Mx My Mz = split cells into Mx by My x Mz (def = 2 2 2)
#                        Mz set to 1 for 2d
#          -x maxlevel = max # of splitting levels (def = 0)
#                        0 = no limit (determined by d)
#                        1 = initial coarse grid, 2 = additional level, etc
#          -d delta = refine grid cells to this size ratio
#                     relative to surf element size (def = 1)
#                     d = 0, split if cell contains a surf, up to maxlevel
#                     d = 1, grid cells will be same size
#                            as any surf element they contain
#                     d = 0.5, grid cells will be half the size
#                              of any surf element they contain
#                     d = 2, grid cells will be 2x the size
#                            of any surf element they contain
#          -f Fx Fy Fz = vector of flow direction (def = 0 0 0)
#            if specified, only surfs with normal against flow are refined
#          -o outfile = SPARTA grid file (def = gdata.tmp)

# NOTE: could enforce that cells without surfs will be no more
#         than 1 level different than neighbors
#       this would propagate refinement away from surf

from __future__ import print_function

import os,sys
from math import sqrt
path = os.environ["SPARTA_PYTHON_TOOLS"]
sys.path.append(path)
from sdata import sdata
from cut2d import cliptest as clip2d
from cut3d import clip as clip3d

# error message

def error(str):
  if str: print("ERROR:",str)
  else: print("Syntax: grid_refine.py -s surffile " + \
        "-b xlo xhi ylo yhi zlo zhi -n Nx Ny Nz -m Mx My Mz " + \
        "-x maxlevel -d delta -f Fx Fy Fz -o outfile")
  sys.exit()

# check a list of surfs for intersection with box
# slist = indices of surfs that intersected with parent cell of this box
# return new sub-list of elements that intersect box
  
def intersect(box,slist):
  newlist = []
  if dim == 2:
    for i in slist:
      p = pts[lines[i][0]]
      q = pts[lines[i][1]]
      cell = (0,box[0],box[1],box[2],box[3])
      if clip2d(p,q,cell): newlist.append(i)
  if dim == 3:
    for i in slist:
      p0 = pts[tris[i][0]]
      p1 = pts[tris[i][1]]
      p2 = pts[tris[i][2]]
      cell = (0,box[0],box[1],box[2],box[3],box[4],box[5])
      if clip3d(cell,p0,p1,p2): newlist.append(i)
  return newlist

# convert ID list into SPARTA ID string
# root is special case
  
def id2str(id):
  if len(id) == 1: return "0"
  string = ""
  for i,value in enumerate(id[1:]):
    string += str(value+1)
    if i < len(id[1:])-1: string += "-"
  return string

# distance between 2 pts

def dist2d(p,q):
  dx = p[0] - q[0]
  dy = p[1] - q[1]
  return sqrt(dx*dx + dy*dy)

def dist3d(p,q):
  dx = p[0] - q[0]
  dy = p[1] - q[1]
  dz = p[2] - q[2]
  return sqrt(dx*dx + dy*dy + dz*dz)

# dot products

def dot2d(a,b):
  dot = a[0]*b[0] + a[1]*b[1]
  return dot

def dot3d(a,b):
  dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
  return dot

# line and tri normal dirs (unnormalized)

def norm2d(a,b):
  x = a[1] - b[1]
  y = b[0] - a[0]
  return (x,y,0)

def norm3d(a,b,c):
  ba = (b[0]-a[0],b[1]-a[1],b[2]-a[2])
  ca = (c[0]-a[0],c[1]-a[1],c[2]-a[2])
  x = ba[1]*ca[2] - ba[2]*ca[1]
  y = ba[2]*ca[0] - ba[0]*ca[2]
  z = ba[0]*ca[1] - ba[1]*ca[0]
  return (x,y,z)

# ---------------------------------------------------------------------
# main program

# read args

surffile = ""
boxxlo = boxylo = boxzlo = 0.0
boxxhi = boxyhi = boxzhi = 1.0
nx = ny = nz = 2
mx = my = mz = 2
maxlevel = 0
delta = 1.0
flowflag = 0
fx = fy = fz = 0.0
outfile = "gdata.tmp"

arg = sys.argv
narg = len(sys.argv)

iarg = 1
while iarg < narg:
  if arg[iarg] == "-s":
    if iarg+2 > narg: error("")
    surffile = arg[iarg+1]
    iarg += 2
  elif arg[iarg] == "-b":
    if iarg+7 > narg: error("")
    boxxlo = float(arg[iarg+1])
    boxxhi = float(arg[iarg+2])
    boxylo = float(arg[iarg+3])
    boxyhi = float(arg[iarg+4])
    boxzlo = float(arg[iarg+5])
    boxzhi = float(arg[iarg+6])
    iarg += 7
  elif arg[iarg] == "-n":
    if iarg+4 > narg: error("")
    nx = int(arg[iarg+1])
    ny = int(arg[iarg+2])
    nz = int(arg[iarg+3])
    iarg += 4
  elif arg[iarg] == "-m":
    if iarg+4 > narg: error("")
    mx = int(arg[iarg+1])
    my = int(arg[iarg+2])
    mz = int(arg[iarg+3])
    iarg += 4
  elif arg[iarg] == "-x":
    if iarg+2 > narg: error("")
    maxlevel = int(arg[iarg+1])
    iarg += 2
  elif arg[iarg] == "-d":
    if iarg+2 > narg: error("")
    delta = float(arg[iarg+1])
    iarg += 2
  elif arg[iarg] == "-f":
    if iarg+4 > narg: error("")
    fx = float(arg[iarg+1])
    fy = float(arg[iarg+2])
    fz = float(arg[iarg+3])
    if fx != 0.0 or fy != 0.0 or fz != 0.0: flowflag = 1
    iarg += 4
  elif arg[iarg] == "-o":
    if iarg+2 > narg: error("")
    outfile = arg[iarg+1]
    iarg += 2
  else: error("")
  
# error check

if not surffile: error("No surface file specified")
if delta < 0.0: error("")
if delta == 0 and maxlevel == 0: error("Delta = 0.0 requires maxlevel > 0")
  
# read surf file via sdata
  
s = sdata(0,surffile)
dim = s.dim
pts = s.surfs[0].points
lines = s.surfs[0].lines
tris = s.surfs[0].triangles

if dim == 2: nz = mz = 1
  
# more error checks

if boxxlo >= boxxhi or boxylo >= boxyhi: error("Bad box size")
if dim == 3 and boxzlo >= boxzhi: error("Bad box size")
if dim == 2 and fz != 0.0: error("Bad Fz value for 2d")
  
# compute size of each surf element
# for tri, size = minimum edge length

sizes = []

if dim == 2:
  for line in lines:
    sizes.append(dist2d(pts[line[0]],pts[line[1]]))
if dim == 3:
  for tri in tris:
    d1 = dist3d(pts[tri[0]],pts[tri[1]])
    d2 = dist3d(pts[tri[1]],pts[tri[2]])
    d3 = dist3d(pts[tri[2]],pts[tri[0]])
    size = min(d1,d2)
    size = min(size,d3)
    sizes.append(size)

# compute surf norms if flowflag set

if flowflag:
  norms = []
  if dim == 2:
    for line in lines:
      norms.append(norm2d(pts[line[0]],pts[line[1]]))
  if dim == 3:
    for tri in tris:
      norms.append(norm3d(pts[tri[0]],pts[tri[1]],pts[tri[2]]))

# initial slist = all surfs
# unless flowflag is set, then exclude those with normals in flow direction

if dim == 2:
  if flowflag:
    slist = []
    for i,norm in enumerate(norms):
      if dot2d((fx,fy),norms[i]) <= 0.0: slist.append(i)
  else: slist = list(range(len(lines)))
if dim == 3:
  if flowflag:
    slist = []
    for i,norm in enumerate(norms):
      if dot3d((fx,fy,fz),norms[i]) <= 0.0: slist.append(i)
  else: slist = list(range(len(tris)))

# generate grid, one level at a time
# queue = list of parent cells to process, remove from front, add to end
#   each entry: id (as list), bbox, surf-list of intersections, nny, nny, nnz
#   compute new intersections from parent intersections
# plist = SPARTA grid file format for saved parent cells
# save no info for child cells, just count them

bbox = (boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi)
queue = [([0],bbox,slist,nx,ny,nz)]
plist = [(0,nx,ny,nz)]
nchild = 0

while len(queue):
  parent = queue.pop(0)
  id = parent[0]
  bbox = parent[1]
  slist = parent[2]
  nnx = parent[3]
  nny = parent[4]
  nnz = parent[5]
  
  index = 0
  for k in range(nnz):
    for j in range(nny):
      for i in range(nnx):
        xlo = i*(bbox[1]-bbox[0])/nnx + bbox[0]
        xhi = (i+1)*(bbox[1]-bbox[0])/nnx + bbox[0]
        if i+1 == nnx: xhi = bbox[1]
        ylo = j*(bbox[3]-bbox[2])/nny + bbox[2]
        yhi = (j+1)*(bbox[3]-bbox[2])/nny + bbox[2]
        if j+1 == nny: yhi = bbox[3]
        zlo = k*(bbox[5]-bbox[4])/nnz + bbox[4]
        zhi = (k+1)*(bbox[5]-bbox[4])/nnz + bbox[4]
        if k+1 == nnz: zhi = bbox[5]
        
        newbox = (xlo,xhi,ylo,yhi,zlo,zhi)
        slistnew = intersect(newbox,slist)

        flag = 0
        if slistnew:
          if delta == 0.0:
            if len(id) < maxlevel: flag = 1
          else:
            cellsize = min(xhi-xlo,yhi-ylo)
            if dim == 3: cellsize = min(cellsize,zhi-zlo)
            for m in slistnew:
              if cellsize/sizes[m] > delta: flag = 1 

        if flag:
          newid = id + [index]
          if dim == 2:
            queue.append((newid,newbox,slistnew,mx,my,1))
            plist.append((id2str(newid),mx,my,1))
          if dim == 3:
            queue.append((newid,newbox,slistnew,mx,my,mz))
            plist.append((id2str(newid),mx,my,mz))
        else: nchild += 1
        
        index += 1

# write out SPARTA grid file
        
print("writing grid file %s ..." % outfile)
fp = open(outfile,"w")

print("# SPARTA grid file produced by refine.py tool", file=fp)
print(file=fp)
print("%d parents" % len(plist), file=fp)
print(file=fp)
print("Parents", file=fp)
print(file=fp)

count = 0
for parent in plist:
  count += 1
  print(count,parent[0],parent[1],parent[2],parent[3], file=fp)
  
fp.close()

# final stats

print(len(plist),"parent cells")
print(nchild,"child cells")
