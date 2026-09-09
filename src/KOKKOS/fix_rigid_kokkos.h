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

#ifdef FIX_CLASS

FixStyle(rigid/kk,FixRigidKokkos)

#else

#ifndef SPARTA_FIX_RIGID_KOKKOS_H
#define SPARTA_FIX_RIGID_KOKKOS_H

#include "fix_rigid.h"

namespace SPARTA_NS {

class FixRigidKokkos : public FixRigid {
 public:
  FixRigidKokkos(class SPARTA *, int, char **);
  ~FixRigidKokkos() {}
  void init();
  void setup();
  void start_of_step();
  void end_of_step();
  void grid_changed();

 private:
  int last_body();
  void host_begin();
  void host_end();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Fix rigid/kk requires compute surf/kk

The compute surf used by fix rigid must be the KOKKOS version, so that
the force/torque tallies are computed by the KOKKOS particle mover.

*/
