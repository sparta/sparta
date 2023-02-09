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

FixStyle(ave/histo/weight/kk,FixAveHistoWeightKokkos)

#else

#ifndef SPARTA_FIX_AVE_HISTO_WEIGHT_KOKKOS_H
#define SPARTA_FIX_AVE_HISTO_WEIGHT_KOKKOS_H

#include "fix_ave_histo_kokkos.h"

namespace SPARTA_NS {

struct TagFixAveHistoWeight_BinVector {};
struct TagFixAveHistoWeight_BinParticles1 {};
struct TagFixAveHistoWeight_BinParticles2 {};
struct TagFixAveHistoWeight_BinParticles3 {};
struct TagFixAveHistoWeight_BinParticles4 {};
struct TagFixAveHistoWeight_BinGridCells1 {};
struct TagFixAveHistoWeight_BinGridCells2 {};
struct TagFixAveHistoWeight_BinParticlesX1 {};
struct TagFixAveHistoWeight_BinParticlesX2 {};
struct TagFixAveHistoWeight_BinParticlesX3 {};
struct TagFixAveHistoWeight_BinParticlesX4 {};
struct TagFixAveHistoWeight_BinParticlesV1 {};
struct TagFixAveHistoWeight_BinParticlesV2 {};
struct TagFixAveHistoWeight_BinParticlesV3 {};
struct TagFixAveHistoWeight_BinParticlesV4 {};

class FixAveHistoWeightKokkos : public FixAveHistoKokkos {
 public:
  FixAveHistoWeightKokkos(class SPARTA *, int, char **);
  virtual ~FixAveHistoWeightKokkos();

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinVector, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinParticles1, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinParticles2, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinParticles3, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinParticles4, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinGridCells1, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinGridCells2, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinParticlesX1, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinParticlesX2, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinParticlesX3, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinParticlesX4, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinParticlesV1, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinParticlesV2, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinParticlesV3, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHistoWeight_BinParticlesV4, const int, minmax_type::value_type&) const;


 private:
  int stridewt;

  using FixAveHisto::bin_one;
  using FixAveHisto::bin_vector;
  using FixAveHisto::bin_particles;
  using FixAveHisto::bin_grid_cells;

  using FixAveHistoKokkos::bin_one;

  // override these methods from ave/histo to use weights
  void bin_scalar(typename minmax_type::value_type&, double);

  void bin_vector(minmax_type&, int, double *, int);
  void bin_particles(minmax_type&, int, int);
  void bin_particles(minmax_type&, double *, int);
  void bin_grid_cells(minmax_type&, DAT::t_float_1d_strided);

  void calculate_weights();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Fix ave/histo/weight value and weight vector lengths do not match

Self-explanatory.

E: Invalid timestep reset for fix ave/histo

Resetting the timestep has invalidated the sequence of timesteps this
fix needs to process.

E: Error writing out histogram data

Something in the output to the file triggered an error.

*/
