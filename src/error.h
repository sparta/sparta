/* ----------------------------------------------------------------------
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

#ifndef SPARTA_ERROR_H
#define SPARTA_ERROR_H

#include "pointers.h"

#include <string>

namespace SPARTA_NS {

class Error : protected Pointers {
 public:
  enum { ERROR_NONE = 0, ERROR_NORMAL = 1, ERROR_ABORT = 2 };

  Error(class SPARTA *);

  void universe_all(const char *, int, const char *);
  void universe_one(const char *, int, const char *);

  void all(const char *, int, const char *);
  void one(const char *, int, const char *);
  void warning(const char *, int, const char *, int = 1);
  void message(const char *, int, const char *, int = 1);
  void done();

  // last error message handling for the library interface

  void set_last_error(const char *msg, int type = ERROR_NORMAL);
  const char *get_last_error() const { return last_error_message.c_str(); }
  int get_last_error_type() const { return last_error_type; }

 private:
  std::string last_error_message;
  int last_error_type = ERROR_NONE;
};

}

#endif
