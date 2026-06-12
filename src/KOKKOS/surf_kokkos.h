/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com
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
  ~SurfKokkos() override;
  void clear_explicit() override;
  void wrap_kokkos();
  void grow(int) override;
  void grow_own(int) override;
  void sync(ExecutionSpace, unsigned int);
  void modify(ExecutionSpace, unsigned int);

  int add_custom(char *, int, int) override;
  void allocate_custom(int) override;
  void reallocate_custom() override;
  void remove_custom(int) override;
  void spread_custom(int) override;
  void spread_inverse_custom(int) override;
  int pack_custom(int, char *) override;
  int unpack_custom(char *, double *) override;

  tdual_line_1d k_lines;
  tdual_tri_1d k_tris;

  tdual_line_1d k_mylines;
  tdual_tri_1d k_mytris;

  DAT::tdual_int_1d k_ewhich,k_eicol,k_edcol;

  tdual_struct_tdual_int_1d_1d k_eivec;
  tdual_struct_tdual_float_1d_1d k_edvec;
  tdual_struct_tdual_int_2d_1d k_eiarray;
  tdual_struct_tdual_float_2d_1d k_edarray;

  tdual_struct_tdual_int_1d_1d k_eivec_local;
  tdual_struct_tdual_float_1d_1d k_edvec_local;
  tdual_struct_tdual_int_2d_1d k_eiarray_local;
  tdual_struct_tdual_float_2d_1d k_edarray_local;
};

}

#endif

/* ERROR/WARNING messages:

*/
