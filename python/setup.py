#!/usr/local/bin/python

"""
setup.py file for DSMC with system MPI library
"""

from distutils.core import setup, Extension

import os, glob
path = os.path.dirname(os.getcwd())

# list of src files for DSMC

libfiles = glob.glob("%s/src/*.cpp" % path)

dsmc_library = Extension("_dsmc",
                         sources = libfiles,
                         define_macros = [("MPICH_IGNORE_CXX_SEEK",1),
                                          ("DSMC_GZIP",1),],
                         # src files for DSMC
                         include_dirs = ["../src"],
                         # additional libs for MPICH on Linux
                         libraries = ["mpich","rt"],
                         # where to find the MPICH lib on Linux
                         library_dirs = ["/usr/local/lib"],
                         # additional libs for MPI on Mac
                         # libraries = ["mpi"],
                         )

setup(name = "dsmc",
      version = "23Sep11",
      author = "Steve Plimpton",
      author_email = "sjplimp@sandia.gov",
      url = "www.sandia.gov/~sjplimp/dsmc.html",
      description = """DSMC library - parallel""",
      py_modules = ["dsmc"],
      ext_modules = [dsmc_library]
      )
