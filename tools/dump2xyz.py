#!/usr/bin/env python

# Script:  dump2xyz.py
# Purpose: convert a SPARTA particle dump file to XYZ format
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# Syntax:  dump2xyz.py dumpfile Nid Ntype Nx Ny Nz xyzfile
#          dumpfile = SPARTA particle dump file
#          Nid,Ntype,Nx,Ny,Nz = columns #s for ID,type,x,y,z
#                               (usually 1,2,3,4,5)
#          xyzfile = new XYZ file

import sys,os
path = os.environ["SPARTA_PYTHON_TOOLS"]
sys.path.append(path)
from dump import dump
from xyz import xyz

if len(sys.argv) != 8:
  raise Exception("Syntax: dump2xyz.py dumpfile Nid Ntype Nx Ny Nz xyzfile")

dumpfile = sys.argv[1]
nid = int(sys.argv[2])
ntype = int(sys.argv[3])
nx = int(sys.argv[4])
ny = int(sys.argv[5])
nz = int(sys.argv[6])
xyzfile = sys.argv[7]

d = dump(dumpfile)
d.map(nid,"id",ntype,"type",nx,"x",ny,"y",nz,"z")
x = xyz(d)
x.one(xyzfile)
