#!/bin/bin/env python

# Script:  surf2stl.py
# Purpose: convert a SPARTA surface file into an STL file (ASCII text)
# Author:  Andrew Hong (Sandia), ayhong at sandia.gov
# Syntax:  surf2stl.py stlfile surffile
#          surffile = write this SPARTA surface file
#          stlfile = read this stereolithography (STL) file
#                    in ASCII (text) format (not binary)

from __future__ import print_function

# error message

def error(str=None):
  if str: print("ERROR:",str)
  else: print("Syntax: surf2stl.py prefix Nfiles step-per-file")
  sys.exit()

def tri_normal(p1,p2,p3):
  N = np.cross(p2-p1, p3-p1);
  N = N / np.linalg.norm(N)
  return N

def convert(fname,outname):

  # get header and find number of points and triangles

  with open(fname) as f:
    head = [next(f) for _ in range(4)];

  npoints = int(head[2].split()[0]);
  ntris = int(head[3].split()[0]);

  # get and store points and trianglulartions

  surftxt = open(fname,"r");

  points = np.zeros((npoints+1,3))
  tris = np.zeros((ntris+1,3));

  lines = surftxt.readlines();

  start = 7;
  for i in range(start,start+npoints):
    line = lines[i].split()
    ind = int(line[0])
    points[ind][0] = float(line[1]);
    points[ind][1] = float(line[2]);
    points[ind][2] = float(line[3]);

  start = npoints+7+3;

  # first triangle may not be 1

  ind = 1
  for i in range(start,start+ntris):
    line = lines[i].split();
    tris[ind][0] = int(line[1]);
    tris[ind][1] = int(line[2]);
    tris[ind][2] = int(line[3]);
    ind = ind + 1;

  # write STL file

  fp = open(outname,"w")

  # header

  print("solid STL",file=fp);


  for i in range(1,ntris+1):
    tri = tris[i];
    p0 = points[int(tri[0])];
    p1 = points[int(tri[1])];
    p2 = points[int(tri[2])];
    n = tri_normal(p0,p1,p2);

    print("  facet normal {0:.6e} {1:.6e} {2:.6e}".format(n[0],n[1],n[2]),file=fp);
    print("    outer loop",file=fp);
    print("      vertex   {0:.6e}  {1:.6e}  {2:.6e}".format(p0[0],p0[1],p0[2]),\
          file=fp);
    print("      vertex   {0:.6e}  {1:.6e}  {2:.6e}".format(p1[0],p1[1],p1[2]),\
          file=fp);
    print("      vertex   {0:.6e}  {1:.6e}  {2:.6e}".format(p2[0],p2[1],p2[2]),\
          file=fp);
    print("    endloop",file=fp);
    print("  endfacet",file=fp);

  print("endsolid vcg",file=fp);
  return

# ----------------------------------------------------------------------
# main program

import sys,re
import numpy as np

inname = sys.argv[1]
outname = sys.argv[2]

ioutname = outname + ".stl"
convert(inname,ioutname)

sys.exit()
