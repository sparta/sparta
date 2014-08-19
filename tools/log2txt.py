#!/usr/bin/env python

# Script:  log2txt.py
# Purpose: extract stats info from SPARTA log file
#          create a text file of numbers in columns, suitable for plotting
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# Syntax:  log2txt.py log.sparta data.txt X Y ...
#          log.sparta = SPARTA log file
#          data.txt = text file to create
#          X Y ... = columns to include (optional), X,Y are stats keywords
#                    if no columns listed, all columns are included

import sys,os
path = os.environ["SPARTA_PYTHON_TOOLS"]
sys.path.append(path)
from log import log

if len(sys.argv) < 3:
  raise StandardError, "Syntax: log2txt.py log.sparta data.txt X Y ..."

logfile = sys.argv[1]
datafile = sys.argv[2]
columns = sys.argv[3:]

lg = log(logfile)
if columns == []:
  lg.write(datafile)
else:
  str = "lg.write(datafile,"
  for word in columns: str += '"' + word + '",'
  str = str[:-1] + ')'
  eval(str)
