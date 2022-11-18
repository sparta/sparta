#!usr/bin/env python

# Script:  implicit_grid.py
# Purpose: create a binary file suitable for reading via the read_isurf command
#          the file is random values on a 2d or 3d grid to
#          infer implicit surfaces from, e.g. a porous material
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# Syntax:  implicit_grid.py dim Nx Ny Nz seed file
#          dim = 2 or 3 for 2d/3d
#          Nx,Ny,Nz = size of corner point grid (for 2d, Nz = 1 is required)
#          seed = random number seed
#          file = name of binary file to ouutput

# Notes:
# generate random grid values from 0 to 255 inclusive on 2d or 3d grid
# explicitly set boundary values to 0 on all 4 or 6 edges/faces
# write to binary file

import sys,random,struct
import numpy as np

def error():
  print("Syntax: implicit_grid.py dim Nx Ny Nz seed file")
  sys.exit()

# command-line args

args = sys.argv[1:]
if len(args) != 6: error()

dim = int(args[0])
nx = int(args[1])
ny = int(args[2])
nz = int(args[3])
seed = int(args[4])
file = args[5]

if dim < 2 or dim > 3: error()
if nx < 3 or ny < 3: error()
if dim == 2 and nz != 1: error()
if dim == 3 and nz < 3: error

# a = 2d or 3d grid with C-style ordering
# to match SPARTA, x varies fastest, y middle, z slowest

if dim == 2: a = np.empty([ny,nx],np.uint8)
else: a = np.empty([nz,ny,ny],np.uint8)

# populate grid
# interior = random ints from 0 to 255 inclusive
# edges/faces = 0

random.seed(seed)

if dim == 2:
  for j in range(ny):
    for i in range(nx):
      if i == 0 or i == nx-1 or j == 0 or j == ny-1: a[j][i] = 0
      else: a[j][i] = random.randint(0,255)
else:
  for k in range(nz):
    for j in range(ny):
      for i in range(nx):
        if i == 0 or i == nx-1 or j == 0 or j == ny-1 or \
           k == 0 or k == nz-1: a[k][j][i] = 0
        else: a[k][j][i] = random.randint(0,255)

# output as binary file
# xyz = 2 or 3 ints, bytes = contents of array in C order
# better to use tobytes() instead of tostring() ??

if dim == 2: xyz = struct.pack("2i",nx,ny)
else: xyz = struct.pack("3i",nx,ny,nz)
bytes = a.tostring()
open(file,'w').write(xyz+bytes)
