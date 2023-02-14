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

#ifdef FIX_CLASS

FixStyle(adapt/kk,FixAdaptKokkos)

#else

#ifndef SPARTA_FIX_ADAPT_KOKKOS_H
#define SPARTA_FIX_ADAPT_KOKKOS_H

#include "fix_adapt.h"

namespace SPARTA_NS {

class FixAdaptKokkos : public FixAdapt {
 public:
  FixAdaptKokkos(class SPARTA *, int, char **);
  ~FixAdaptKokkos() {}
  void end_of_step();
};

}

#endif
#endif
