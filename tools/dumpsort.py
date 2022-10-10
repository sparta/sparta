#!/usr/bin/env python

# Script:  dumpsort.py
# Purpose: sort the snapshots in a SPARTA particle dump file by particle ID
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# Syntax:  dumpsort.py oldfile N newfile
#          oldfile = old SPARTA dump file
#          N = column # for atom ID (usually 1)
#          newfile = new sorted SPARTA dump file

import sys,os
path = os.environ["SPARTA_PYTHON_TOOLS"]
sys.path.append(path)
from dump import dump

if len(sys.argv) != 4:
  raise Exception("Syntax: dumpsort.py oldfile N newfile")

oldfile = sys.argv[1]
ncolumn = int(sys.argv[2])
newfile = sys.argv[3]

d = dump(oldfile)
d.map(ncolumn,"id")
d.sort()
d.write(newfile)
