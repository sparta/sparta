# ----------------------------------------------------------------------
#   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
#   http://sparta.sandia.gov
#   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
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

import sys,traceback,types
from ctypes import *

class sparta:
  def __init__(self,name="",cmdargs=None):

    # load libsparta.so by default
    # if name = "g++", load libsparta_g++.so

    try:
      if not name: self.lib = CDLL("libsparta.so",RTLD_GLOBAL)
      else: self.lib = CDLL("libsparta_%s.so" % name,RTLD_GLOBAL)
    except:
      type,value,tb = sys.exc_info()
      traceback.print_exception(type,value,tb)
      raise Exception("Could not load SPARTA dynamic library")

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
      # self.spa = self.lib.sparta_open_no_mpi(0,None)

  def __del__(self):
    if self.spa: self.lib.sparta_close(self.spa)

  def close(self):
    self.lib.sparta_close(self.spa)
    self.spa = None

  def file(self,file):
    file = file.encode('utf-8')
    self.lib.sparta_file(self.spa,file)

  def command(self,cmd):
    cmd = cmd.encode('utf-8')
    self.lib.sparta_command(self.spa,cmd)

  def extract_global(self,name,type):
    name = name.encode('utf-8')
    if type == 0:
      self.lib.sparta_extract_global.restype = POINTER(c_int)
    elif type == 1:
      self.lib.sparta_extract_global.restype = POINTER(c_double)
    else: return None
    ptr = self.lib.sparta_extract_global(self.spa,name)
    return ptr[0]

  def extract_compute(self,id,style,type):
    style = style.encode('utf-8')
    if type == 0:
      self.lib.sparta_extract_compute.restype = POINTER(c_double)
      ptr = self.lib.sparta_extract_compute(self.spa,id,style,type)
      return ptr[0]
    if type == 1:
      self.lib.sparta_extract_compute.restype = POINTER(c_double)
      ptr = self.lib.sparta_extract_compute(self.spa,id,style,type)
      return ptr
    if type == 2:
      self.lib.sparta_extract_compute.restype = POINTER(POINTER(c_double))
      ptr = self.lib.sparta_extract_compute(self.spa,id,style,type)
      return ptr
    return None

  # free memory for 1 double or 1 vector of doubles via sparta_free()
  # for vector, must copy nlocal returned values to local c_double vector
  # memory was allocated by library interface function

  def extract_variable(self,name,type):
    name = name.encode('utf-8')
    if type == 0:
      self.lib.sparta_extract_variable.restype = POINTER(c_double)
      ptr = self.lib.sparta_extract_variable(self.spa,name)
      result = ptr[0]
      self.lib.sparta_free(ptr)
      return result
    if type == 1:
      self.lib.sparta_extract_global.restype = POINTER(c_int)
      nlocalptr = self.lib.sparta_extract_global(self.spa,b"nplocal")
      nlocal = nlocalptr[0]
      result = (c_double*nlocal)()
      self.lib.sparta_extract_variable.restype = POINTER(c_double)
      ptr = self.lib.sparta_extract_variable(self.spa,name)
      for i in range(nlocal): result[i] = ptr[i]
      self.lib.sparta_free(ptr)
      return result
    return None
