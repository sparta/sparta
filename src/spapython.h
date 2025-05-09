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

#ifndef SPARTA_PYTHON_H
#define SPARTA_PYTHON_H

#include "pointers.h"

namespace SPARTA_NS {

class PythonInterface {
 public:
  virtual ~PythonInterface() noexcept(false) {}
  virtual void command(int, char **) = 0;
  virtual void invoke_function(int, char *) = 0;
  virtual int find(const char *) = 0;
  virtual int variable_match(const char *, const char *, int) = 0;
  virtual char *long_string(int ifunc) = 0;
  virtual int execute_string(char *) = 0;
  virtual int execute_file(char *) = 0;
  virtual bool has_minimum_version(int major, int minor) = 0;
};

class Python : protected Pointers {
 public:
  Python(class SPARTA *);
  ~Python() override;

  void command(int, char **);
  void invoke_function(int, char *);
  int find(const char *);
  int variable_match(const char *, const char *, int);
  char *long_string(int ifunc);
  int execute_string(char *);
  int execute_file(char *);
  bool has_minimum_version(int major, int minor);

  bool is_enabled() const;
  void init();
  static void finalize();

 private:
  PythonInterface *impl;
};

}

#endif
