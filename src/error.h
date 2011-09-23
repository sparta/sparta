/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#ifndef DSMC_ERROR_H
#define DSMC_ERROR_H

#include "pointers.h"

namespace DSMC_NS {

class Error : protected Pointers {
 public:
  Error(class DSMC *);

  void universe_all(const char *, int, const char *);
  void universe_one(const char *, int, const char *);

  void all(const char *, int, const char *);
  void one(const char *, int, const char *);
  void warning(const char *, int, const char *, int = 1);
  void message(const char *, int, char *, int = 1);
  void done();
};

}

#endif
