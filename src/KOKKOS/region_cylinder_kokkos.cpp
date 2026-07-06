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

#include "region_cylinder_kokkos.h"
#include "sparta_masks.h"
#include "memory_kokkos.h"
#include "particle_kokkos.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

RegCylinderKokkos::RegCylinderKokkos(SPARTA *sparta, int narg, char **arg) :
  RegCylinder(sparta, narg, arg)
{
  kokkos_flag = 1;
}

/* ---------------------------------------------------------------------- */

RegCylinderKokkos::RegCylinderKokkos(SPARTA *sparta) :
  RegCylinder(sparta)
{
  copy = 1;
}

/* ---------------------------------------------------------------------- */

RegCylinderKokkos::~RegCylinderKokkos()
{
}

/* ---------------------------------------------------------------------- */

void RegCylinderKokkos::match_all_kokkos(DAT::tdual_int_1d k_match_in)
{
  d_match = k_match_in.view_device();
  ParticleKokkos* particleKK = (ParticleKokkos*) particle;
  particleKK->sync(Device, PARTICLE_MASK);
  d_particles = particleKK->k_particles.view_device();
  int nlocal = particle->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<TagRegCylinderMatchAll>(0,nlocal),*this);
  copymode = 0;
  k_match_in.modify_device();
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void RegCylinderKokkos::operator()(TagRegCylinderMatchAll, const int &i) const
{
  const double x_tmp = d_particles[i].x[0];
  const double y_tmp = d_particles[i].x[1];
  const double z_tmp = d_particles[i].x[2];
  d_match[i] = match_kokkos(x_tmp,y_tmp,z_tmp);
}
