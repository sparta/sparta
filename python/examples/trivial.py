#!/usr/bin/env python -i
# preceeding line should have path for Python on your machine

# trivial.py
# Purpose: run a SPARTA input script via Python
# Syntax:  trivial.py in.sparta
#          in.sparta = any SPARTA input script

import sys

# parse command line

argv = sys.argv
if len(argv) != 2:
  print "Syntax: trivial.py in.sparta"
  sys.exit()

infile = sys.argv[1]

me = 0
# uncomment if running in parallel via Pypar
#import pypar
#me = pypar.rank()
#nprocs = pypar.size()

from sparta import sparta
spa = sparta()

# run infile all at once

spa.file(infile)

# run infile one line at a time

#lines = open(infile,'r').readlines()
#for line in lines: spa.command(line)

# uncomment if running in parallel via Pypar
#print "Proc %d out of %d procs has" % (me,nprocs), lmp
#pypar.finalize()
