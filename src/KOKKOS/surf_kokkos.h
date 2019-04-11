/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov
   Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_SURF_KOKKOS_H
#define SPARTA_SURF_KOKKOS_H

#include "stdio.h"
#include "surf.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

class SurfKokkos : public Surf {
 public:
  SurfKokkos(class SPARTA *);
  ~SurfKokkos();
  void wrap_kokkos();
  void sync(ExecutionSpace, unsigned int);
  void modify(ExecutionSpace, unsigned int);
  void grow();
  void grow_own();

  tdual_line_1d k_lines;
  tdual_tri_1d k_tris;

  tdual_line_1d k_mylines;
  tdual_tri_1d k_mytris;

 private:

};

}

#endif

/* ERROR/WARNING messages:

*/
