# ----------------------------------------------------------------------
#   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
#   http://sparta.sandia.gov
#   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
#   Sandia National Laboratories
#
#   Copyright (2012) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under 
#   the GNU General Public License.
#
#   See the README file in the top-level SPARTA directory.
# -------------------------------------------------------------------------

# Python wrapper on SPARTA library via ctypes

import types,sys
from ctypes import *

class sparta:
  def __init__(self,name="",cmdargs=None):

    # load libsparta.so by default
    # if name = "g++", load libsparta_g++.so
    
    try:
      if not name: self.lib = CDLL("libsparta.so")
      else: self.lib = CDLL("libsparta_%s.so" % name)
    except:
      type,value,tb = sys.exc_info()
      traceback.print_exception(type,value,tb)
      raise OSError,"Could not load SPARTA dynamic library"

    # create an instance of SPARTA
    # don't know how to pass an MPI communicator from PyPar
    # no_mpi call lets SPARTA use MPI_COMM_WORLD
    # cargs = array of C strings from args

    if cmdargs:
      cmdargs.insert(0,"sparta.py")
      narg = len(cmdargs)
      cargs = (c_char_p*narg)(*cmdargs)
      self.spa = c_void_p()
      self.lib.sparta_open_no_mpi(narg,cargs,byref(self.spa))
    else:
      self.spa = c_void_p()
      self.lib.sparta_open_no_mpi(0,None,byref(self.spa))
      # could use just this if SPARTA lib interface supported it
      # self.lmp = self.lib.lammps_open_no_mpi(0,None)

  def __del__(self):
    if self.spa: self.lib.sparta_close(self.spa)

  def close(self):
    self.lib.sparta_close(self.spa)
    self.spa = None

  def file(self,file):
    self.lib.sparta_file(self.spa,file)

  def command(self,cmd):
    self.lib.sparta_command(self.spa,cmd)
