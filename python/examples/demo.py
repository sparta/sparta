#!/usr/bin/env python -i
# preceeding line should have path for Python on your machine

# demo.py
# Purpose: illustrate use of many library interface commands
# Syntax:  demo.py
#          uses in.demo as SPARTA input script

import sys

# parse command line

argv = sys.argv
if len(argv) != 1:
  print("Syntax: demo.py")
  sys.exit()

me = 0
# uncomment if running in parallel via Pypar
#import pypar
#me = pypar.rank()
#nprocs = pypar.size()

from sparta import sparta

spa = sparta()

# test out various library functions after running in.demo

spa.file(b"in.demo")

if me == 0: print("\nPython output:")

nparticles = spa.extract_global(b"nplocal",0)
dt = spa.extract_global(b"dt",1)
fnum = spa.extract_global(b"fnum",1)
print("Nparticles, dt, fnum =",nparticles,dt,fnum)

temp = spa.extract_compute(b"temp",0,0)
print("Temperature from compute =",temp)

lx = spa.extract_variable(b"lx",0)
print("Box-length lx from equal-style variable =",lx)

xc = spa.extract_variable(b"xc",1)
vy = spa.extract_variable(b"vy",1)
print("X coord from particle-style variable =",xc[0],xc[nparticles-1])
print("Vy component coord from particle-style variable =",vy[0],vy[nparticles-1])

# uncomment if running in parallel via Pypar
#print "Proc %d out of %d procs has" % (me,nprocs), spa
#pypar.finalize()
