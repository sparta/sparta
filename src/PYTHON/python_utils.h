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

#ifndef SPARTA_PYTHON_UTILS_H
#define SPARTA_PYTHON_UTILS_H

#include <Python.h>

namespace SPARTA_NS {

namespace PyUtils {

  class GIL {
    PyGILState_STATE gstate;

   public:
    GIL() : gstate(PyGILState_Ensure()) {}
    ~GIL() { PyGILState_Release(gstate); }
  };

  static void Print_Errors()
  {
    PyErr_Print();
    PyErr_Clear();
  }

}

}

#endif
