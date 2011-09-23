# ----------------------------------------------------------------------
#   DSMC - Sandia parallel DSMC code
#   www.sandia.gov/~sjplimp/dsmc.html
#   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
#   Sandia National Laboratories
#
#   Copyright (2011) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under 
#   the GNU General Public License.
#
#   See the README file in the top-level DSMC directory.
# -------------------------------------------------------------------------

# Python wrapper on DSMC library via ctypes

import types
from ctypes import *

DSMCINT = 0
DSMCDOUBLE = 1
DSMCIPTR = 2
DSMCDPTR = 3
DSMCDPTRPTR = 4

class dsmc:
  def __init__(self,args=None):

    # attempt to load parallel library first, serial library next
    # could provide caller a flag to choose which library to load
    
    try:
      self.lib = CDLL("_dsmc.so")
    except:
      try:
        self.lib = CDLL("_dsmc_serial.so")
      except:
        raise StandardError,"Could not load DSMC dynamic library"

    # create an instance of DSMC
    # don't know how to pass an MPI communicator from PyPar
    # no_mpi call lets DSMC use MPI_COMM_WORLD
    # cargs = array of C strings from args
    
    if args:
      args.insert(0,"dsmc.py")
      narg = len(args)
      cargs = (c_char_p*narg)(*args)
      self.dsmc = c_void_p()
      self.lib.dsmc_open_no_mpi(narg,cargs,byref(self.dsmc))
    else:
      self.dsmc = c_void_p()
      self.lib.dsmc_open_no_mpi(0,None,byref(self.dsmc))
      # could use just this if DSMC lib interface supported it
      # self.dsmc = self.lib.dsmc_open_no_mpi(0,None)

  def __del__(self):
    if self.dsmc: self.lib.dsmc_close(self.dsmc)

  def close(self):
    self.lib.dsmc_close(self.dsmc)
    self.dsmc = None

  def file(self,file):
    self.lib.dsmc_file(self.dsmc,file)

  def command(self,cmd):
    self.lib.dsmc_command(self.dsmc,cmd)
