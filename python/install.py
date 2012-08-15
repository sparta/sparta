#!/usr/bin/env python

instructions = """copy SPARTA shared library src/libsparta.so and sparta.py to system dirs
Syntax: python install.py [libdir] [pydir]
    libdir = target dir for src/libsparta.so, default = /usr/local/lib, or the first
        item in LD_LIBRARY_PATH if it doesn't exist.
    pydir = target dir for sparta.py, default = Python site-packages, via distutils."""

import sys, shutil, os

if len(sys.argv) > 3:
  print instructions
  sys.exit()

# verify that our user-specified path is in LD_LIBRARY_PATH
# since if not, the install won't work
  
libdir = "/usr/local/lib"
if "LD_LIBRARY_PATH" in os.environ:
  libpaths = os.environ['LD_LIBRARY_PATH'].split(':')
  if not libdir in libpaths: libdir = libpaths[0]

pydir = False
try:
  libdir = sys.argv[1]
  pydir = sys.argv[2]
except IndexError:
  pass

# copy the C library into place

shutil.copy('../src/libsparta.so', libdir)

# if user-specified, copy sparta.py into directory
# else invoke setup from Distutils to add to site-packages

if pydir:
  shutil.copy('../python/sparta.py', pydir)
  sys.exit()

from distutils.core import setup

os.chdir('../python')

setup(name = "sparta",
      version = "15Aug12",
      author = "Steve Plimpton",
      author_email = "sjplimp@sandia.gov",
      url = "http://sparta.sandia.gov",
      description = """SPARTA DSMC library""",
      py_modules = ["sparta"])
