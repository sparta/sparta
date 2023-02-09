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

CommandStyle(read_surf/kk,ReadSurfKokkos)

#else

#ifndef SPARTA_READ_SURF_KOKKOS_H
#define SPARTA_READ_SURF_KOKKOS_H

#include "read_surf.h"

namespace SPARTA_NS {

class ReadSurfKokkos : public ReadSurf {
 public:
  ReadSurfKokkos(class SPARTA *);
  ~ReadSurfKokkos() {}
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
