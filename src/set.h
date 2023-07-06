/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(set,Set)

#else

#ifndef SPARTA_SET_H
#define SPARTA_SET_H

#include "pointers.h"

namespace SPARTA_NS {

class Set : protected Pointers {

 public:
  Set(class SPARTA *);
  void command(int, char **);

 private:
  int mode,groupbit;
  int ctype,csize,cindex,ccol;
  int vindex,vstyle;
  char *cname,*vname;
  
  class Mixture *mixture;
  class Region *region;
  class Variable *variable;

  int set_particle(double, double *);
  int set_grid(double, double *);
  int set_surf(double, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
