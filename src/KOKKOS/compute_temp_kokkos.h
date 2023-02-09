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

#ifdef COMPUTE_CLASS

ComputeStyle(temp/kk,ComputeTempKokkos)

#else

#ifndef SPARTA_COMPUTE_TEMP_KOKKOS_H
#define SPARTA_COMPUTE_TEMP_KOKKOS_H

#include "compute_temp.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

class ComputeTempKokkos : public ComputeTemp {
 public:
  ComputeTempKokkos(class SPARTA *, int, char **);
  virtual ~ComputeTempKokkos() {}
  double compute_scalar();

  KOKKOS_INLINE_FUNCTION
  void operator()(const int&, double&) const;

 private:
  t_particle_1d d_particles;
  t_species_1d d_species;

  double compute_scalar_kokkos();
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
