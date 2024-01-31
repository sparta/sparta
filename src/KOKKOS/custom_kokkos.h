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

CommandStyle(custom/k,CustomKokkos)

#else

#ifndef SPARTA_CUSTOM_KOKKOS_H
#define SPARTA_CUSTOM_KOKKOS_H

#include "custom.h"

namespace SPARTA_NS {

class CustomKokkos : public Custom {

 public:
  CustomKokkos(class SPARTA *);
  void command(int, char **);

 private:
  int set_particle(double, double *);
  int set_grid(double, double *);
  int set_surf(double, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
