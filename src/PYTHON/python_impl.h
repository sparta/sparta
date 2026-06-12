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

#ifndef SPARTA_PYTHON_IMPL_H
#define SPARTA_PYTHON_IMPL_H

#include "spapython.h"
#include "pointers.h"

namespace SPARTA_NS {

class PythonImpl : protected Pointers, public PythonInterface {

 public:
  PythonImpl(class SPARTA *);
  ~PythonImpl() override;
  void command(int, char **) override;
  void invoke_function(int, char *, double *) override;
  int find(const char *) override;
  int function_match(const char *, const char *, int) override;
  int wrapper_match(const char *, const char *, int, int *) override;
  char *long_string(int) override;
  int execute_string(char *) override;
  int execute_file(char *) override;
  bool has_minimum_version(int major, int minor) override;
  static void finalize();

 private:
  void *pyMain;

  struct PyFunc {
    char *name;
    int ninput, noutput;
    int *itype, *ivarflag;
    int *ivalue;
    double *dvalue;
    char **svalue;
    int *internal_var;        // stores per-arg index of internal variable
    int otype;
    char *ovarname;
    char *longstr;
    int length_longstr;
    void *pFunc;
  };

  PyFunc *pfuncs;
  int nfunc;

  int create_entry(char *, int, int, int, char **, char *, char *);
  void deallocate(int);
};

}

#endif
