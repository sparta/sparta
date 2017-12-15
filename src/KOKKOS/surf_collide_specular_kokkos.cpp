/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "surf_collide_specular_kokkos.h"
#include "surf.h"
#include "surf_react.h"
#include "modify.h"
#include "math_extra.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

SurfCollideSpecularKokkos::SurfCollideSpecularKokkos(SPARTA *sparta, int narg, char **arg) :
  SurfCollideSpecular(sparta, narg, arg)
{
  k_nsingle = DAT::tdual_int_scalar("SurfCollide:nsingle");
  d_nsingle = k_nsingle.view<DeviceType>();
  h_nsingle = k_nsingle.h_view;
}

SurfCollideSpecularKokkos::SurfCollideSpecularKokkos(SPARTA *sparta) :
  SurfCollideSpecular(sparta)
{
  // ID and style
  // ID must be all alphanumeric chars or underscores

  int narg = 2;
  const char* arg[] = {"sc_kk_specular_copy","specular"};

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,"Surf_collide ID must be alphanumeric or "
                 "underscore characters");

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  vector_flag = 1;
  size_vector = 2;
    
  nsingle = ntotal = 0;

  copy = 0;

  if (narg != 2) error->all(FLERR,"Illegal surf_collide specular command");

  k_nsingle = DAT::tdual_int_scalar("SurfCollide:nsingle");
  d_nsingle = k_nsingle.view<DeviceType>();
  h_nsingle = k_nsingle.h_view;
}
