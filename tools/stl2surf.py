#!/bin/bin/env python

# Script:  stl2surf.py
# Purpose: convert an STL file (ASCII text) into a SPARTA surface file
#          and warn if surface is not watertight
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# Syntax:  stl2surf.py stlfile surffile
#          stlfile = read this stereolithography (STL) file
#                    in ASCII (text) format (not binary)
#          surffile = write this SPARTA surface file

# NOTE: process vertices in text format so no precision or round-off issues

from __future__ import print_function

# error message

def error(str=None):
  if str: print("ERROR:",str)
  else: print("Syntax: stl2surf.py stlfile surffile")
  sys.exit()

# ----------------------------------------------------------------------
# main program

import sys,re

if len(sys.argv) != 3: error()

stlfile = sys.argv[1]
surffile = sys.argv[2]

# parse STL file into triangles and triangle vertices
# tritxt = list of text between facet and endfacet
# triverts = list of 3 vertices per triangle, in text format

stltxt = open(stlfile,"r").read()

match = re.search("^.*\n",stltxt)
if not match: error("STL file %s has incorrect format" % stlfile)
words = match.group().split()
if words[0] != "solid":
  error("STL file %s has incorrect format" % stlfile)
if len(words) >= 2: name = words[1]
else: name = ""

tritxt = re.split("endfacet\s+facet",stltxt)
print("# of triangles in STL file:",len(tritxt))

pattern = ".*vertex\s+(\S+)\s+(\S+)\s+(\S+)"
triverts = []

for one in tritxt:
  vertices = re.findall(pattern,one)
  if len(vertices) != 3:
    print("Triangle record:")
    print(one)
    error("Invalid triangle vertices")
  triverts.append(vertices)

# build list of unique vertices via hash
# unique: key = vertex 3-tuple, value = index in verts
# verts = list of unique vertices, each as 3-tuple
# tris = list of 3 vertex indices for each triangle

unique = {}
verts = []
tris = []

for vert3 in triverts:
  if vert3[0] in unique: v0 = unique[vert3[0]]
  else:
    v0 = len(verts)
    verts.append(vert3[0])
    unique[vert3[0]] = v0
  if vert3[1] in unique: v1 = unique[vert3[1]]
  else:
    v1 = len(verts)
    verts.append(vert3[1])
    unique[vert3[1]] = v1
  if vert3[2] in unique: v2 = unique[vert3[2]]
  else:
    v2 = len(verts)
    verts.append(vert3[2])
    unique[vert3[2]] = v2
  tris.append((v0,v1,v2))

# print SPARTA surface file

fp = open(surffile,"w")

if name:
  print("# SPARTA surface file, from STL file %s with name %s\n" % \
      (stlfile,name),file=fp)
else:
  print("# SPARTA surface file, from STL file\n",stlfile,file=fp)

print(len(verts),"points",file=fp)
print(len(tris),"triangles",file=fp)

print("\nPoints\n",file=fp)
for i,vert in enumerate(verts):
  print(i+1,vert[0],vert[1],vert[2],file=fp)

print("\nTriangles\n",file=fp)
for i,tri in enumerate(tris):
  print(i+1,tri[0]+1,tri[1]+1,tri[2]+1,file=fp)

fp.close()
  
# stats to screen

print("# of vertices in SPARTA file:",len(verts))
print("# of triangles in SPARTA file:",len(tris))

# warn if not a watertight object
# watertight = each edge is part of exactly 2 triangles
#   and edge direction is different in 2 triangles

ehash = {}
dup = 0
unmatch = 0

for vert3 in triverts:
  edge = (vert3[0],vert3[1])
  if edge in ehash:
    dup += 1
    dupedge = edge
  else: ehash[edge] = 1
  edge = (vert3[1],vert3[2])
  if edge in ehash:
    dup += 1
    dupedge = edge
  else: ehash[edge] = 1
  edge = (vert3[2],vert3[0])
  if edge in ehash:
    dup += 1
    dupedge = edge
  else: ehash[edge] = 1

for edge in ehash:
  invedge = (edge[1],edge[0])
  if invedge not in ehash:
    unmatch += 1
    unmatchedge = edge
    
if dup or unmatch:
  print("WARNING: surface is not watertight")
  print("Duplicate edge count:",dup)
  print("Unmatched edge count:",unmatch)
  print("One duplicate edge:",dupedge)
  print("One unmatched edge:",unmatchedge)
