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

#include "math.h"
#include "surf_collide_specular.h"
#include "surf.h"
#include "surf_react.h"
#include "modify.h"
#include "math_extra.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

SurfCollideSpecular::SurfCollideSpecular(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg < 2 || narg > 3) error->all(FLERR,"Illegal surf_collide specular command");

  // optional args

  noslip_flag = 0;

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"noslip") == 0) {
      noslip_flag = 1;
      iarg += 1;
    } else
      error->all(FLERR,"Illegal surf_collide specular command");
  }
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   ip = particle with current x = collision pt, current v = incident v
   isurf = index of surface element
   norm = surface normal unit vector
   isr = index of reaction model if >= 0, -1 for no chemistry
   ip = reset to NULL if destroyed by chemistry
   return jp = new particle if created by chemistry
   return reaction = index of reaction (1 to N) that took place, 0 = no reaction
   resets particle(s) to post-collision outward velocity
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideSpecular::
collide(Particle::OnePart *&ip, double &,
        int isurf, double *norm, int isr, int &reaction)
{
  nsingle++;

  // if surface chemistry defined, attempt reaction
  // reaction = 1 to N for which reaction took place, 0 for none
  // velreset = 1 if reaction reset post-collision velocity, else 0

  Particle::OnePart iorig;
  Particle::OnePart *jp = NULL;
  reaction = 0;
  int velreset = 0;

  if (isr >= 0) {
    if (modify->n_surf_react) memcpy(&iorig,ip,sizeof(Particle::OnePart));
    reaction = surf->sr[isr]->react(ip,isurf,norm,jp,velreset);
    if (reaction) surf->nreact_one++;
  }

  // specular or noslip reflection for each particle
  // only if SurfReact did not already reset velocities
  // also both particles need to trigger any fixes
  //   to update per-particle properties which depend on
  //   temperature of the particle, e.g. fix vibmode and fix ambipolar
  // NOTE: not doing this for this specular model,
  //   since temperature does not change, would need to add a twall arg

  if (ip) {
    if (!velreset) {
      if (noslip_flag) {

        // noslip reflection
        // reflect incident v, all three components

        MathExtra::negate3(ip->v);
      } else {

        // specular reflection

        MathExtra::reflect3(ip->v,norm);
      }
    }

    //if (modify->n_update_custom) {
    //  int i = ip - particle->particles;
    //  modify->update_custom(i,twall,twall,twall,vstream);
    //}
  }

  if (jp) {
    if (!velreset) {
      if (noslip_flag) {

        // noslip reflection
        // reflect incident v, all three components.

        MathExtra::negate3(jp->v);
      } else {

        // specular reflection

        MathExtra::reflect3(jp->v,norm);
      }
    }

    //if (modify->n_update_custom) {
    //  int j = jp - particle->particles;
    //  modify->update_custom(j,twall,twall,twall,vstream);
    //}
  }

  // call any fixes with a surf_react() method
  // they may reset j to -1, e.g. fix ambipolar
  //   in which case newly created j is deleted

  if (reaction && modify->n_surf_react) {
    int i = -1;
    if (ip) i = ip - particle->particles;
    int j = -1;
    if (jp) j = jp - particle->particles;
    modify->surf_react(&iorig,i,j);
    if (jp && j < 0) {
      jp = NULL;
      particle->nlocal--;
    }
  }

  return jp;
}

/* ----------------------------------------------------------------------
   wrapper on specular() method to perform collision for a single particle
   pass in 1 flag (optional) to match command-line args for style specular
   flags, coeffs can be NULL
   called by SurfReactAdsorb
------------------------------------------------------------------------- */

void SurfCollideSpecular::wrapper(Particle::OnePart *p, double *norm,
                                  int *flags, double *coeffs)
{
  if (flags)
    noslip_flag = flags[0];

  MathExtra::reflect3(p->v,norm);
}

/* ----------------------------------------------------------------------
   return flags and coeffs for this SurfCollide instance to caller
------------------------------------------------------------------------- */

void SurfCollideSpecular::flags_and_coeffs(int *flags, double *coeffs)
{
  flags[0] = noslip_flag;
}
