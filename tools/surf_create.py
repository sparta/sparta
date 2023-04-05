#!/usr/bin/env python

# Script:  surf_create.py
# Purpose: create a simple triangulated (3d) or line-segmented (2d) surfaces
#          output file can contain one or more surface objects
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# Syntax:  surf_create.py sfile style args style args ...
#          sfile = surface file to create
#          style = sphere or box or spikysphere (3d objects)
#            sphere args = x y z rad N
#              x,y,z = center point
#              rad = radius
#              N = NxN squares per face of cube embedded in sphere
#            box args = xlo ylo zlo xhi yhi zhi Nx Ny Nz
#              xlo,ylo,zlo = lower left corner point
#              xhi,yhi,zhi = upper right corner point
#              Nx,Ny,Nz = # of divisions in each dimension of box
#            spikysphere args = x y z rmin rmax N
#              x,y,z = center point
#              rmin,rmax = lo/hi bounds of random radius of each sphere surf pt
#              N = NxN squares per face of cube embedded in sphere
#          style = circle or rect or tri or spikycircle (2d objects)
#            circle args = x y rad N
#              x,y = center point
#              rad = radius
#              N = # of line segments around circumference of circle
#            rect args = xlo ylo xhi yhi Nx Ny
#              xlo,ylo = lower left corner point
#              xhi,yhi = upper right corner point
#              Nx,Ny = # of line segments in each dimension
#            tri args = x1 y1 x2 y2 x3 y3 N1 N2 N3
#              x1,y1 = first point of triangle
#              x2,y2 = second point of triangle
#              x3,y3 = third point of triangle
#              N1,N2,N3 = # of line segments for each side of triangle
#            spikycircle args = x y rmin rmax N
#              x,y = center point
#              rmin,rmax = lo/hi bounds of random radius of each circle surf pt
#              N = # of line segements around circumference of circle

# NOTE: cannot mix 2d and 3d objects

import sys,os
path = os.environ["SPARTA_PYTHON_TOOLS"]
sys.path.append(path)
from sdata import sdata

# error message

def error():
  print("Syntax: surf_create.py sfile style argsx style args ...")
  sys.exit()

# ---------------------------------------------------------------------
# main program

args = sys.argv
narg = len(args)

if narg < 3: error()

sfile = sys.argv[1]

s = sdata()
id = "a"

iarg = 2
while iarg < narg:
  if args[iarg] == "sphere":
    if iarg+6 > narg: error()
    s.sphere(id,float(args[iarg+1]),float(args[iarg+2]),
             float(args[iarg+3]),float(args[iarg+4]),int(args[iarg+5]))
    iarg += 6
  elif args[iarg] == "box":
    if iarg+10 > narg: error()
    s.box(id,float(args[iarg+1]),float(args[iarg+2]),
          float(args[iarg+3]),float(args[iarg+4]),
          float(args[iarg+5]),float(args[iarg+6]),
          int(args[iarg+7]),int(args[iarg+8]),int(args[iarg+9]))
    iarg += 10
  elif args[iarg] == "spikysphere":
    if iarg+7 > narg: error()
    s.spikysphere(id,float(args[iarg+1]),float(args[iarg+2]),
                  float(args[iarg+3]),float(args[iarg+4]),
                  float(args[iarg+5]),int(args[iarg+6]))
    iarg += 7
  elif args[iarg] == "circle":
    if iarg+5 > narg: error()
    s.circle(id,float(args[iarg+1]),float(args[iarg+2]),
             float(args[iarg+3]),int(args[iarg+4]))
    iarg += 5
  elif args[iarg] == "rect":
    if iarg+7 > narg: error()
    s.rect(id,float(args[iarg+1]),float(args[iarg+2]),
           float(args[iarg+3]),float(args[iarg+4]),
           int(args[iarg+5]),int(args[iarg+6]))
    iarg += 7
  elif args[iarg] == "tri":
    if iarg+10 > narg: error()
    s.tri(id,float(args[iarg+1]),float(args[iarg+2]),
          float(args[iarg+3]),float(args[iarg+4]),
          float(args[iarg+5]),float(args[iarg+6]),
          int(args[iarg+7]),int(args[iarg+8]),int(args[iarg+9]))
    iarg += 10
  elif args[iarg] == "spikycircle":
    if iarg+6 > narg: error()
    s.spikycircle(ud,float(args[iarg+1]),float(args[iarg+2]),
                  float(args[iarg+3]),float(args[iarg+4]),nt(args[iarg+5]))
    iarg += 6
  else: error()
  id += "a"
  
s.write(sfile)
