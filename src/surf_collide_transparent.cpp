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
#include "string.h"
#include "surf_collide_transparent.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

SurfCollideTransparent::
SurfCollideTransparent(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg != 2) error->all(FLERR,"Illegal surf_collide transparent command");

  transparent = 1;
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   ip = particle with current x = collision pt, current v = incident v
   norm = surface normal unit vector
   isr = index of reaction model if >= 0, -1 for no chemistry
   simply return ip = NULL to delete particle
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideTransparent::
collide(Particle::OnePart *&ip, double *, double &, int, int &)
{
  nsingle++;

  return NULL;
}
