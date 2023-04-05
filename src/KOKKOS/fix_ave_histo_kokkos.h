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

FixStyle(ave/histo/kk, FixAveHistoKokkos)

#else

#ifndef SPARTA_FIX_AVE_HISTO_KOKKOS_H
#define SPARTA_FIX_AVE_HISTO_KOKKOS_H

#include "fix_ave_histo.h"
#include "kokkos_type.h"
#include "grid_kokkos.h"

namespace SPARTA_NS
{

struct TagFixAveHisto_BinVector {};
struct TagFixAveHisto_BinParticles1 {};
struct TagFixAveHisto_BinParticles2 {};
struct TagFixAveHisto_BinParticles3 {};
struct TagFixAveHisto_BinParticles4 {};
struct TagFixAveHisto_BinGridCells1 {};
struct TagFixAveHisto_BinGridCells2 {};
struct TagFixAveHisto_BinParticlesX1 {};
struct TagFixAveHisto_BinParticlesX2 {};
struct TagFixAveHisto_BinParticlesX3 {};
struct TagFixAveHisto_BinParticlesX4 {};
struct TagFixAveHisto_BinParticlesV1 {};
struct TagFixAveHisto_BinParticlesV2 {};
struct TagFixAveHisto_BinParticlesV3 {};
struct TagFixAveHisto_BinParticlesV4 {};

namespace FixKokkosDetails {

template<class ValueType, class OutputDeviceType>
inline
Kokkos::View<ValueType*, OutputDeviceType>
mirror_view_from_raw_host_array(ValueType* x, const int size, const int stride)
{
  typedef typename OutputDeviceType::memory_space out_mem_space;
  typedef Kokkos::View<ValueType*, OutputDeviceType> out_view_type;
  typedef Kokkos::View<ValueType*, typename out_view_type::array_layout,
                       Kokkos::HostSpace> host_view_type;
  host_view_type h_x("h_x", size);
  for (int i=0, j=0; i<size; i++, j+=stride) h_x(i) = x[j];
  out_view_type d_x = Kokkos::create_mirror_view(out_mem_space(), h_x);
  Kokkos::deep_copy(d_x, h_x);
  return d_x;
}

}

class FixAveHistoKokkos : public FixAveHisto
{
public:

  typedef Kokkos::MinMax<double,SPAHostType> minmax_type;
  typedef minmax_type::value_type mm_value_type;

  FixAveHistoKokkos(class SPARTA *, int, char **);
  virtual ~FixAveHistoKokkos();
  void init();
  void setup();
  void end_of_step();
  double compute_vector(int);
  double compute_array(int, int);

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinVector, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinParticles1, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinParticles2, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinParticles3, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinParticles4, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinGridCells1, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinGridCells2, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinParticlesX1, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinParticlesX2, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinParticlesX3, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinParticlesX4, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinParticlesV1, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinParticlesV2, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinParticlesV3, const int, minmax_type::value_type&) const;

  KOKKOS_INLINE_FUNCTION void
  operator()(TagFixAveHisto_BinParticlesV4, const int, minmax_type::value_type&) const;

protected:

  DAT::tdual_float_1d k_stats;
  DAT::t_float_1d d_stats;

  DAT::tdual_float_1d k_bin;
  DAT::t_float_1d d_bin;

  t_particle_1d d_particles;
  DAT::t_int_2d d_s2g;

  int index;
  int stride;
  DAT::t_float_1d_strided d_values;
  GridKokkos* grid_kk;

  // data used by ave/histo/weight/kk
  DAT::t_float_1d_strided d_weights;

  // methods
  using FixAveHisto::bin_one;
  using FixAveHisto::bin_vector;
  using FixAveHisto::bin_particles;
  using FixAveHisto::bin_grid_cells;

  virtual void bin_scalar(mm_value_type&, double);
  virtual void bin_vector(minmax_type&, int, double *, int);
  virtual void bin_particles(minmax_type&, int, int);
  virtual void bin_particles(minmax_type&, double *, int);
  virtual void bin_grid_cells(minmax_type&, DAT::t_float_1d_strided);

  virtual void calculate_weights() {}

  void options(int, int, char **);
  bigint nextvalid();

  /* ----------------------------------------------------------------------
     bin a single value with weight
     ----------------------------------------------------------------------- */
  KOKKOS_INLINE_FUNCTION
  void
  bin_one(mm_value_type& mm_v, double value, double weight) const
  {
    auto nbins = d_bin.extent(0);
    if (value < mm_v.min_val) mm_v.min_val = value;
    if (value > mm_v.max_val) mm_v.max_val = value;
    if (value < lo) {
      if (beyond == 0 /*IGNORE*/) {
        Kokkos::atomic_add(&d_stats(1), weight);
        return;
      } else {
        Kokkos::atomic_add(&d_bin(0), weight);
      }
    } else if (value > hi) {
      if (beyond == 0 /*IGNORE*/) {
        Kokkos::atomic_add(&d_stats(1), weight);
        return;
      } else {
        Kokkos::atomic_add(&d_bin(nbins-1), weight);
      }
    } else {
      int ibin = static_cast<int>((value - lo) * bininv);
      ibin = MIN(ibin, nbins-1);
      if (beyond == 2 /*EXTRA*/) ibin++;
      Kokkos::atomic_add(&d_bin(ibin), weight);
    }
    Kokkos::atomic_add(&d_stats(0), weight);
  }

  KOKKOS_INLINE_FUNCTION
  void
  bin_one(mm_value_type& mm_v, double value) const
  {
    bin_one(mm_v, value, 1.);
  }

};
} // namespace SPARTA_NS

#endif
#endif

    /* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute ID for fix ave/histo does not exist

Self-explanatory.

E: Fix ID for fix ave/histo does not exist

Self-explanatory.

E: Fix ave/histo input is invalid compute

Self-explanatory.

E: Fix ave/histo input is invalid fix

Self-explanatory.

E: Fix ave/histo input is invalid variable

Self-explanatory.

E: Fix ave/histo inputs are not all global, peratom, or local

All inputs in a single fix ave/histo command must be of the
same style.

E: Fix ave/histo cannot input per-atom values in scalar mode

Self-explanatory.

E: Fix ave/histo cannot input local values in scalar mode

Self-explanatory.

E: Fix ave/histo compute does not calculate a global scalar

Self-explanatory.

E: Fix ave/histo compute does not calculate a global vector

Self-explanatory.

E: Fix ave/histo compute vector is accessed out-of-range

Self-explanatory.

E: Fix ave/histo compute does not calculate a global array

Self-explanatory.

E: Fix ave/histo compute array is accessed out-of-range

Self-explanatory.

E: Fix ave/histo compute does not calculate per-atom values

Self-explanatory.

E: Fix ave/histo compute does not calculate a per-atom vector

Self-explanatory.

E: Fix ave/histo compute does not calculate a per-atom array

Self-explanatory.

E: Fix ave/histo compute does not calculate local values

Self-explanatory.

E: Fix ave/histo compute does not calculate a local vector

Self-explanatory.

E: Fix ave/histo compute does not calculate a local array

Self-explanatory.

E: Fix ave/histo fix does not calculate a global scalar

Self-explanatory.

E: Fix ave/histo fix does not calculate a global vector

Self-explanatory.

E: Fix ave/histo fix vector is accessed out-of-range

Self-explanatory.

E: Fix for fix ave/histo not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/histo is
requesting a value on a non-allowed timestep.

E: Fix ave/histo fix does not calculate a global array

Self-explanatory.

E: Fix ave/histo fix array is accessed out-of-range

Self-explanatory.

E: Fix ave/histo fix does not calculate per-atom values

Self-explanatory.

E: Fix ave/histo fix does not calculate a per-atom vector

Self-explanatory.

E: Fix ave/histo fix does not calculate a per-atom array

Self-explanatory.

E: Fix ave/histo fix does not calculate local values

Self-explanatory.

E: Fix ave/histo fix does not calculate a local vector

Self-explanatory.

E: Fix ave/histo fix does not calculate a local array

Self-explanatory.

E: Variable name for fix ave/histo does not exist

Self-explanatory.

E: Error writing file header

Something in the output to the file triggered an error.

E: Invalid timestep reset for fix ave/histo

Resetting the timestep has invalidated the sequence of timesteps this
fix needs to process.

E: Error writing out histogram data

Something in the output to the file triggered an error.

E: Cannot open fix ave/histo file %s

The specified file cannot be opened.  Check that the path and name are
correct.

*/
