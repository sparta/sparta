#!/usr/bin/env python -i

# Script:  logplot.py
# Purpose: use GnuPlot to plot two columns from a SPARTA log file
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# Syntax:  logplot.py logfile X Y
#          logfile = SPARTA log file
#          X,Y = plot Y versus X where X,Y are stats keywords
#          once plot appears, you are in Python interpreter, type C-D to exit

import sys,os
path = os.environ["SPARTA_PYTHON_TOOLS"]
sys.path.append(path)
from olog import olog
from gnu import gnu

if len(sys.argv) != 4:
  raise Exception("Syntax: logplot.py logfile X Y")

logfile = sys.argv[1]
xlabel = sys.argv[2]
ylabel = sys.argv[3]

lg = olog(logfile)
x,y = lg.get(xlabel,ylabel)
g = gnu()
g.plot(x,y)
print("Type Ctrl-D to exit Python")
