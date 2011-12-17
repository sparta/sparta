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

#ifndef DSMC_DOMAIN_H
#define DSMC_DOMAIN_H

#include "pointers.h"

namespace DSMC_NS {

class Domain : protected Pointers {
 public:
  int box_exist;                    // 0 = not yet created, 1 = exists
  int dimension;                    // 2,3

  double boxlo[3],boxhi[3];         // box global bounds
  double xprd,yprd,zprd;                 // global box dimensions
  double prd[3];                    // array form of dimensions

  Domain(class DSMC *);
  ~Domain();
  void init() {}
  void set_initial_box();
  void set_global_box();
  void print_box(const char *);
};

}

#endif

/* ERROR/WARNING messages:

E: Box bounds are invalid

The box boundaries specified in the read_data file are invalid.  The
lo value must be less than the hi value for all 3 dimensions.

*/
