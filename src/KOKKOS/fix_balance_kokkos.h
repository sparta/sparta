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

FixStyle(balance/kk,FixBalanceKokkos)

#else

#ifndef SPARTA_FIX_BALANCE_KOKKOS_H
#define SPARTA_FIX_BALANCE_KOKKOS_H

#include "fix_balance.h"

namespace SPARTA_NS {

class FixBalanceKokkos : public FixBalance {
 public:
  FixBalanceKokkos(class SPARTA *, int, char **);
  ~FixBalanceKokkos() {}
  void end_of_step();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Cannot use non-rcb fix balance with a grid cutoff

This is because the load-balancing will generate a partitioning
of cells to processors that is dispersed and which will not work
with a grid cutoff >= 0.0.

*/
