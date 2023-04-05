#!usr/bin/env python

# Script:  jagged2d.py
# Purpose: create a 2d jagged surf file to test SPARTA surface options
#          created in x orientation in a box of size 1 with origin at 0
#          read surf command can tranlate/scale as desired
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# Syntax:  jagged2d.py Nspike Nper delta sfile
#          Nspike = number of spikes as seen from -x side (upwind in flow)
#          Nper = number of line segments on each side of each spike
#          delta = extra thickness of 3 line segment closure on +x side
#          sfile = surface data file to output

# Notes:
# add 3 additional lines to create a closed surf between first/last spike

from __future__ import print_function

import sys,random
import numpy as np

# generate nper line segments between p1 and p2

def genlines(p1,p2):
  for i in range(nper):
    x1 = p1[0] + (i+0.0)/nper * (p2[0]-p1[0])
    y1 = p1[1] + (i+0.0)/nper * (p2[1]-p1[1])
    x2 = p1[0] + (i+1.0)/nper * (p2[0]-p1[0])
    y2 = p1[1] + (i+1.0)/nper * (p2[1]-p1[1])
    if i == nper-1:
      x2 = p2[0]; y2 = p2[1]
    pts.append((x1,y1))
    pts.append((x2,y2))
    lines.append((len(pts)-2,len(pts)-1))

# main program

args = sys.argv[1:]
if len(args) != 4:
  print("Syntax: jagged2d.py Nspike Nper delta sfile")
  sys.exit()
  
nspike = int(args[0])
nper = int(args[1])
delta = float(args[2])
datafile = args[3]

# create line segments for all spikes
# traverse jagged outline from one end to the other
#   so that order of pts in each line is correct for flow volume

pts = []
lines = []

for i in range(nspike):
  p1 = (1.0,(i+0.0)/nspike)
  p2 = (0.0,(i+0.5)/nspike)
  genlines(p1,p2)
  
  p1 = (0.0,(i+0.5)/nspike)
  if i == nspike-1: p2 = (1.0,1.0)
  else: p2 = (1.0,(i+1.0)/nspike)
  genlines(p1,p2)

# add 3 line segments to close the surface

pts.append((1.0+delta,0.0))
pts.append((1.0,0.0))
lines.append((len(pts)-2,len(pts)-1))

pts.append((1.0,1.0))
pts.append((1.0+delta,1.0))
lines.append((len(pts)-2,len(pts)-1))

pts.append((1.0+delta,1.0))
pts.append((1.0+delta,0.0))
lines.append((len(pts)-2,len(pts)-1))

# write out surf file

print("Writing %s with %d points, %d lines" % (datafile,len(pts),len(lines)))

fp = open(datafile,'w')

print("surf file from jagged2d.py", file=fp)
print(file=fp)
print(len(pts),"points", file=fp)
print(len(lines),"lines", file=fp)
print(file=fp)
print("Points\n", file=fp)
for i,pt in enumerate(pts):
  print(i+1,pt[0],pt[1], file=fp)
print(file=fp)
print("Lines\n", file=fp)
for i,line in enumerate(lines):
  print(i+1,1,line[0]+1,line[1]+1, file=fp)
    
fp.close()
