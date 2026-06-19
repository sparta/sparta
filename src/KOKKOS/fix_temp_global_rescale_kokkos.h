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

FixStyle(temp/global/rescale/kk,FixTempGlobalRescaleKokkos)

#else

#ifndef SPARTA_FIX_TEMP_GLOBAL_RESCALE_KOKKOS_H
#define SPARTA_FIX_TEMP_GLOBAL_RESCALE_KOKKOS_H

#include "fix_temp_global_rescale.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

struct TagFixTempGlobalRescale_reduce{};
struct TagFixTempGlobalRescale_scale{};

class FixTempGlobalRescaleKokkos : public FixTempGlobalRescale {
 public:
  FixTempGlobalRescaleKokkos(class SPARTA *, int, char **);
  virtual ~FixTempGlobalRescaleKokkos() {}
  void end_of_step() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixTempGlobalRescale_reduce, const int&, double&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixTempGlobalRescale_scale, const int&) const;

 private:
  double vscale;

  t_particle_1d d_particles;
  t_species_1d d_species;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/
