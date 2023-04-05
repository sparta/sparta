#!/usr/bin/env python

# Script:  surf_transform.py
# Purpose: morph a surf file into another surf file
#          with same transformation options as read_surf command
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# Syntax:  surf_transform.py infile outfile keyword args ...
#          infile,outfile = input and output surf files
#          keywords = origin, trans, atrans, scale, rotate, invert
#          origin Ox Oy Oz
#          trans Dx Dy Dz
#          atrans Ax Ay Az
#          scale Sx Sy Sz
#          rotate theta Rx Ry Rz
#          invert (no args)

from __future__ import print_function

import os,sys,math
path = os.environ["SPARTA_PYTHON_TOOLS"]
sys.path.append(path)
from sdata import sdata

# error message

def error(str):
  if str: print("ERROR:",str)
  else: print("Syntax: surf_transform.py infile outfile keyword args ...")
  sys.exit()

# write SPARTA surf file
  
def write_surfs(dim,file):
  fp = open(file,"w")
  print("%s file converted via surf_transform.py\n" % infile, file=fp)
  
  print(len(pts),"points", file=fp)
  if lines: print(len(lines),"lines", file=fp)
  if tris: print(len(tris),"triangles", file=fp)

  print("\nPoints\n", file=fp)
  if dim == 2:
    for i,pt in enumerate(pts):
      print(i+1,pt[0],pt[1], file=fp)
  else:
    for i,pt in enumerate(pts):
      print(i+1,pt[0],pt[1],pt[2], file=fp)

  if lines:
    print("\nLines\n", file=fp)
    for i,line in enumerate(lines):
      print(i+1,line[0]+1,line[1]+1, file=fp)
    
  if tris:
    print("\nTriangles\n", file=fp)
    for i,tri in enumerate(tris):
      print(i+1,tri[0]+1,tri[1]+1,tri[2]+1, file=fp)

  fp.close()

# translate origin and pts by dx,dy,dz
  
def translate(dx,dy,dz):
  origin[0] += dx
  origin[1] += dy
  if dim == 3: origin[2] += dz
  for pt in pts:
    pt[0] += dx
    pt[1] += dy
    if dim == 3: pt[2] += dz

# move origin to ax,ay,ax and translate pts by dx,dy,dz
  
def atranslate(ax,ay,az):
  dx = ax - origin[0]
  dy = ay - origin[1]
  if dim == 3: dz = az - origin[2]
  origin[0] = ax
  origin[1] = ay
  if dim == 3: origin[2] = az
  for pt in pts:
    pt[0] += dx
    pt[1] += dy
    if dim == 3: pt[2] += dz

# scale pts by sx,sy,sz around origin
# for 2d, do not reset x[2] to avoid epsilon change

def scale(dim,sx,sy,sz):
  for pt in pts:
    pt[0] = sx*(pt[0]-origin[0]) + origin[0];
    pt[1] = sy*(pt[1]-origin[1]) + origin[1];
    if dim == 3: pt[2] = sz*(pt[2]-origin[2]) + origin[1];

# rotate pts by theta around rx,ry,rz vector from origin
# for 2d, do not reset x[2] to avoid epsilon change
    
def rotate(dim,theta,rx,ry,rz):
  PI = 3.14159265358979323846
  theta *= PI/180.0
  r = [rx,ry,rz]
  norm3(r)
  q = axisangle_to_quat(r,theta)
  rotmat = quat_to_mat(q)

  d = 3*[0.0]
  for pt in pts:
    d[0] = pt[0] - origin[0];
    d[1] = pt[1] - origin[1];
    d[2] = pt[2] - origin[2];
    dnew = matvec(rotmat,d)
    pt[0] = dnew[0] + origin[0];
    pt[1] = dnew[1] + origin[1];
    if dim == 3: pt[2] = dnew[2] + origin[2];
  
# invert vertex ordering within each line or tri
# this flips direction of surface normal
    
def invert(dim):
  if dim == 2:
    for line in lines:
      tmp = line[1]
      line[0] = line[1]
      line[1] = tmp
  if dim == 3:
    for tri in tris:
      tmp = tri[1]
      tri[1] = tri[2]
      tri[2] = tmp

# normalize a 3 vector

def norm3(x):
  scale = 1.0/math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
  x[0] *= scale
  x[1] *= scale
  x[2] *= scale

def axisangle_to_quat(v,angle):
  halfa = 0.5*angle;
  sina = math.sin(halfa);
  quat = 4*[0.0]
  quat[0] = math.cos(halfa);
  quat[1] = v[0]*sina;
  quat[2] = v[1]*sina;
  quat[3] = v[2]*sina;
  return quat

def quat_to_mat(quat):
  w2 = quat[0]*quat[0];
  i2 = quat[1]*quat[1];
  j2 = quat[2]*quat[2];
  k2 = quat[3]*quat[3];
  twoij = 2.0*quat[1]*quat[2];
  twoik = 2.0*quat[1]*quat[3];
  twojk = 2.0*quat[2]*quat[3];
  twoiw = 2.0*quat[1]*quat[0];
  twojw = 2.0*quat[2]*quat[0];
  twokw = 2.0*quat[3]*quat[0];

  mat = [[0,0,0],[0,0,0],[0,0,0]]
  mat[0][0] = w2+i2-j2-k2;
  mat[0][1] = twoij-twokw;
  mat[0][2] = twojw+twoik;
  mat[1][0] = twoij+twokw;
  mat[1][1] = w2-i2+j2-k2;
  mat[1][2] = twojk-twoiw;
  mat[2][0] = twoik-twojw;
  mat[2][1] = twojk+twoiw;
  mat[2][2] = w2-i2-j2+k2;
  return mat

def matvec(a,x):
  b = 3*[0.0]
  b[0] = a[0][0]*x[0] + a[0][1]*x[1] + a[0][2]*x[2]
  b[1] = a[1][0]*x[0] + a[1][1]*x[1] + a[1][2]*x[2]
  b[2] = a[2][0]*x[0] + a[2][1]*x[1] + a[2][2]*x[2]
  return b

# ---------------------------------------------------------------------
# main program

arg = sys.argv
narg = len(sys.argv)

if narg < 4: error("")
infile = arg[1]
outfile = arg[2]

# read surf file via sdata
  
s = sdata(0,infile)
dim = s.dim
pts = s.surfs[0].points
lines = s.surfs[0].lines
tris = s.surfs[0].triangles

# perform transform operations, one at a time

origin = [0.0,0.0,0.0]

iarg = 3
while iarg < narg:
  if arg[iarg] == "origin":
    if iarg+4 > narg: error("")
    origin[0] = float(arg[iarg+1])
    origin[1] = float(arg[iarg+2])
    origin[2] = float(arg[iarg+3])
    iarg += 4
  elif arg[iarg] == "trans":
    if iarg+4 > narg: error("")
    tx = float(arg[iarg+1])
    ty = float(arg[iarg+2])
    tz = float(arg[iarg+3])
    translate(tx,ty,tz)
    iarg += 4
  elif arg[iarg] == "atrans":
    if iarg+4 > narg: error("")
    ax = float(arg[iarg+1])
    ay = float(arg[iarg+2])
    az = float(arg[iarg+3])
    atranslate(ax,ay,az)
    iarg += 4
  elif arg[iarg] == "scale":
    if iarg+4 > narg: error("")
    sx = float(arg[iarg+1])
    sy = float(arg[iarg+2])
    sz = float(arg[iarg+3])
    scale(dim,sx,sy,sz)
    iarg += 4
  elif arg[iarg] == "rotate":
    if iarg+5 > narg: error("")
    theta = float(arg[iarg+1])
    rx = float(arg[iarg+2])
    ry = float(arg[iarg+3])
    rz = float(arg[iarg+4])
    rotate(dim,theta,rx,ry,rz)
    iarg += 5
  elif arg[iarg] == "invert":
    invert(dim)
    iarg += 1
  else: error("")

write_surfs(dim,outfile)
