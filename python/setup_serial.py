#!/usr/local/bin/python

"""
setup_serial.py file for DSMC with dummy serial MPI library
"""

from distutils.core import setup, Extension

import os, glob
path = os.path.dirname(os.getcwd())

# list of src files for DSMC and MPI STUBS

libfiles = glob.glob("%s/src/*.cpp" % path) + \
           glob.glob("%s/src/STUBS/*.cpp" % path)

dsmc_library = Extension("_lammps_serial",
                         sources = libfiles,
                         define_macros = [("MPICH_IGNORE_CXX_SEEK",1),
                                          ("DSMC_GZIP",1),],
                         # src files for DSMC and MPI STUBS
                         include_dirs = ["../src", "../src/STUBS"]
                         )

setup(name = "dsmc_serial",
      version = "26Oct10",
      author = "Steve Plimpton",
      author_email = "sjplimp@sandia.gov",
      url = "www.sandia.gov/~sjplimp/dsmc.html",
      description = """DSMC library - serial""",
      py_modules = ["dsmc"],
      ext_modules = [dsmc_library]
      )
