#!usr/bin/env python

# Script:  jagged3d.py
# Purpose: create a 3d jagged surf file to test SPARTA surface options
#          created in x orientation in a box of size 1 with origin at 0
#          read surf command can tranlate/scale as desired
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# Syntax:  jagged3d.py Nyspike Nzspike Nxlayer Nyslice Nzslice delta sfile
#          Nyspike,Nzspike = 2d array of spikes as seen from -x side
#                            (upwind in flow)
#          Nxlayer = number of layers in each spike
#          Nyslice,Nzslice = vertical slices in each spike in each dimension
#          delta = extra thickness of box-like triangle closure on +x side
#          sfile = surface data file to output

# Notes:
# add additional tris to create a closed surf between all the spikes

from __future__ import print_function

import sys,random
import numpy as np

# generate small tris withing large (p1,p2,p3) tri
# p1 is point at x = 0, other 2 points are at x = 1
# order of p1,p2,p3 is with flow volume normal
# nslice = # of vertical slices in tri
# nxlayer = # of layers in tri

def gentris(p1,p2,p3,nslice):
  for j in range(nslice):
    
    # p4 = 1st point on line between p2 and p3
    
    x = p2[0] + (j+0.0)/nslice * (p3[0]-p2[0])
    y = p2[1] + (j+0.0)/nslice * (p3[1]-p2[1])
    z = p2[2] + (j+0.0)/nslice * (p3[2]-p2[2])
    p4 = (x,y,z)

    # p5 = 2nd point on line between p2 and p3

    x = p2[0] + (j+1.0)/nslice * (p3[0]-p2[0])
    y = p2[1] + (j+1.0)/nslice * (p3[1]-p2[1])
    z = p2[2] + (j+1.0)/nslice * (p3[2]-p2[2])
    if j == nslice-1: p5 = p3
    else: p5 = (x,y,z)
    
    for i in range(nxlayer):
      
      if i == 0:     # one tri per slice
        
        x1 = p1[0]
        y1 = p1[1]
        z1 = p1[2]
        x2 = p1[0] + (i+1.0)/nxlayer * (p4[0]-p1[0])
        y2 = p1[1] + (i+1.0)/nxlayer * (p4[1]-p1[1])
        z2 = p1[2] + (i+1.0)/nxlayer * (p4[2]-p1[2])
        x3 = p1[0] + (i+1.0)/nxlayer * (p5[0]-p1[0])
        y3 = p1[1] + (i+1.0)/nxlayer * (p5[1]-p1[1])
        z3 = p1[2] + (i+1.0)/nxlayer * (p5[2]-p1[2])
        
        if i == nxlayer-1:
          x2 = p4[0]; y2 = p4[1]; z2 = p4[2]
          x3 = p5[0]; y3 = p5[1]; z3 = p5[2]

        pts.append((x1,y1,z1))
        pts.append((x2,y2,z2))
        pts.append((x3,y3,z3))
        tris.append((len(pts)-3,len(pts)-2,len(pts)-1))
        
      else:     # one quad per slice = 2 tris

        x1 = p1[0] + (i+0.0)/nxlayer * (p4[0]-p1[0])
        y1 = p1[1] + (i+0.0)/nxlayer * (p4[1]-p1[1])
        z1 = p1[2] + (i+0.0)/nxlayer * (p4[2]-p1[2])
        x2 = p1[0] + (i+1.0)/nxlayer * (p4[0]-p1[0])
        y2 = p1[1] + (i+1.0)/nxlayer * (p4[1]-p1[1])
        z2 = p1[2] + (i+1.0)/nxlayer * (p4[2]-p1[2])
        x3 = p1[0] + (i+1.0)/nxlayer * (p5[0]-p1[0])
        y3 = p1[1] + (i+1.0)/nxlayer * (p5[1]-p1[1])
        z3 = p1[2] + (i+1.0)/nxlayer * (p5[2]-p1[2])
        x4 = p1[0] + (i+0.0)/nxlayer * (p5[0]-p1[0])
        y4 = p1[1] + (i+0.0)/nxlayer * (p5[1]-p1[1])
        z4 = p1[2] + (i+0.0)/nxlayer * (p5[2]-p1[2])
        
        if i == nxlayer-1:
          x2 = p4[0]; y2 = p4[1]; z2 = p4[2]
          x3 = p5[0]; y3 = p5[1]; z3 = p5[2]

        pts.append((x1,y1,z1))
        pts.append((x2,y2,z2))
        pts.append((x3,y3,z3))
        pts.append((x4,y4,z4))
        tris.append((len(pts)-4,len(pts)-3,len(pts)-2))  # pts 1,2,3
        tris.append((len(pts)-4,len(pts)-2,len(pts)-1))  # pts 1,3,4
        
# main program
    
args = sys.argv[1:]
if len(args) != 7:
  print("Syntax: jagged3d.py Nyspike Nzspike Nxlayer Nyslice Nzslice delta sfile")
  sys.exit()
  
nyspike = int(args[0])
nzspike = int(args[1])
nxlayer = int(args[2])
nyslice = int(args[3])
nzslice = int(args[4])
delta = float(args[5])
datafile = args[6]

# create triangles for all spikes, 4 large spike triangles per spike
# insure order of pts in each tri is correct for flow volume normal

pts = []
tris = []

for iz in range(nzspike):
  for iy in range(nyspike):

    # 5 pts in spike
    
    p1 = (0.0,(iy+0.5)/nyspike,(iz+0.5)/nzspike)   # tip at x = 0
    p2 = (1.0,(iy+0.0)/nyspike,(iz+0.0)/nzspike)
    p3 = (1.0,(iy+1.0)/nyspike,(iz+0.0)/nzspike)
    p4 = (1.0,(iy+0.0)/nyspike,(iz+1.0)/nzspike)
    p5 = (1.0,(iy+1.0)/nyspike,(iz+1.0)/nzspike)

    # avoid round-off for pts at upper boundaries
    
    if iy == nyspike-1: p3 = (1.0,1.0,(iz+0.0)/nzspike)
    if iz == nzspike-1: p4 = (1.0,(iy+0.0)/nyspike,1.0)
    if iy == nyspike-1 and iz == nzspike-1: p5 = (1.0,1.0,1.0)
    elif iy == nyspike-1: p5 = (1.0,1.0,(iz+1.0)/nzspike)
    elif iz == nzspike-1: p5 = (1.0,(iy+1.0)/nyspike,1.0)

    # 4 large tris per spike
    
    gentris(p1,p2,p4,nzslice)    # lower y triangle
    gentris(p1,p5,p3,nzslice)    # upper y triangle
    gentris(p1,p3,p2,nyslice)    # lower z triangle
    gentris(p1,p4,p5,nyslice)    # upper z triangle

# add triangles to close the surface
# 4 small, skinny faces on yz perimeter, one large cap at x = 1+delta
# Nyspike*Nyslice spacing in y, Nzspike*Nzslice spacing in z
# be careful to avoid round-off differences in points
#   always do double loop over spikes and slices
# order 4 points in each rectangle for flow volume normal
# 2 triangles per rectangle

# y = 0 face and y = 1 face

for y in [0.0,1.0]:
  for i in range(nzspike):
    lo = (i+0.0)/nzspike
    hi = (i+1.0)/nzspike
    if i == nzspike-1: hi = 1.0

    for j in range(nzslice):
      zlo = lo + (j+0.0)/nzslice * (hi-lo)
      zhi = lo + (j+1.0)/nzslice * (hi-lo)
      if j == nzslice-1: zhi = hi

      if y == 0.0:
        p1 = (1.0,y,zlo)
        p2 = (1.0+delta,y,zlo)
        p3 = (1.0+delta,y,zhi)
        p4 = (1.0,y,zhi)
      else:
        p1 = (1.0,y,zlo)
        p2 = (1.0,y,zhi)
        p3 = (1.0+delta,y,zhi)
        p4 = (1.0+delta,y,zlo)

      pts.append(p1); pts.append(p2); pts.append(p3); pts.append(p4)
      tris.append((len(pts)-4,len(pts)-3,len(pts)-2))  # pts 1,2,3
      tris.append((len(pts)-4,len(pts)-2,len(pts)-1))  # pts 1,3,4
    
# z = 0 face and z = 1 face

for z in [0.0,1.0]:
  for i in range(nyspike):
    lo = (i+0.0)/nyspike
    hi = (i+1.0)/nyspike
    if i == nyspike-1: hi = 1.0

    for j in range(nyslice):
      ylo = lo + (j+0.0)/nyslice * (hi-lo)
      yhi = lo + (j+1.0)/nyslice * (hi-lo)
      if j == nyslice-1: yhi = hi

      if z == 1.0:
        p1 = (1.0,ylo,z)
        p2 = (1.0+delta,ylo,z)
        p3 = (1.0+delta,yhi,z)
        p4 = (1.0,yhi,z)
      else:
        p1 = (1.0,ylo,z)
        p2 = (1.0,yhi,z)
        p3 = (1.0+delta,yhi,z)
        p4 = (1.0+delta,ylo,z)

      pts.append(p1); pts.append(p2); pts.append(p3); pts.append(p4)
      tris.append((len(pts)-4,len(pts)-3,len(pts)-2))  # pts 1,2,3
      tris.append((len(pts)-4,len(pts)-2,len(pts)-1))  # pts 1,3,4

# x = 1+delta cap

for iz in range(nzspike):
  zlo = (iz+0.0)/nzspike
  zhi = (iz+1.0)/nzspike
  if iz == nzspike-1: zhi = 1.0

  for izz in range(nzslice):
    zzlo = zlo + (izz+0.0)/nzslice * (zhi-zlo)
    zzhi = zlo + (izz+1.0)/nzslice * (zhi-zlo)
    if izz == nzslice-1: zzhi = zhi

    for iy in range(nyspike):
      ylo = (iy+0.0)/nyspike
      yhi = (iy+1.0)/nyspike
      if iy == nyspike-1: yhi = 1.0

      for iyy in range(nyslice):
        yylo = ylo + (iyy+0.0)/nyslice * (yhi-ylo)
        yyhi = ylo + (iyy+1.0)/nyslice * (yhi-ylo)
        if iyy == nyslice-1: yyhi = yhi

        p1 = (1.0+delta,yylo,zzlo)
        p2 = (1.0+delta,yyhi,zzlo)
        p3 = (1.0+delta,yyhi,zzhi)
        p4 = (1.0+delta,yylo,zzhi)

        pts.append(p1); pts.append(p2); pts.append(p3); pts.append(p4)
        tris.append((len(pts)-4,len(pts)-3,len(pts)-2))  # pts 1,2,3
        tris.append((len(pts)-4,len(pts)-2,len(pts)-1))  # pts 1,3,4
    
# write out surf file

print("Writing %s with %d points, %d tris" % (datafile,len(pts),len(tris)))

fp = open(datafile,'w')

print("surf file from jagged3d.py", file=fp)
print(file=fp)
print(len(pts),"points", file=fp)
print(len(tris),"triangles", file=fp)
print(file=fp)
print("Points\n", file=fp)
for i,pt in enumerate(pts):
  print(i+1,pt[0],pt[1],pt[2], file=fp)
print(file=fp)
print("Triangles\n", file=fp)
for i,tri in enumerate(tris):
  print(i+1,1,tri[0]+1,tri[1]+1,tri[2]+1, file=fp)
    
fp.close()
