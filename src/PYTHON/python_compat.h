/* -*- c++ -*- ----------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_PYTHON_COMPAT_H
#define SPARTA_PYTHON_COMPAT_H

#include <Python.h>

#if PY_VERSION_HEX < 0x030600f0
#error Python version 3.6 or later is required by SPARTA
#endif

// Wrap API differences between platforms using macros
#if defined(_MSC_VER) || defined(__MINGW32__)
#define PY_INT_FROM_LONG(X) PyLong_FromLongLong(X)
#define PY_INT_AS_LONG(X) PyLong_AsLongLong(X)
#define PY_LONG_FROM_STRING(X) std::stoll(X)
#else
#define PY_INT_FROM_LONG(X) PyLong_FromLong(X)
#define PY_INT_AS_LONG(X) PyLong_AsLong(X)
#define PY_LONG_FROM_STRING(X) std::stol(X)
#endif

#endif
