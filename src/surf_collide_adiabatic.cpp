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
#include "surf_collide_adiabatic.h"
#include "surf.h"
#include "surf_react.h"
#include "update.h"
#include "particle.h"
#include "modify.h"
#include "comm.h"
#include "math_const.h"
#include "math_extra.h"
#include "error.h"
#include "random_mars.h"
#include "random_knuth.h"

using namespace SPARTA_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

SurfCollideAdiabatic::SurfCollideAdiabatic(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg != 2) error->all(FLERR,"Illegal surf_collide adiabatic command");

  // initialize RNG for isotropic reflection

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfCollideAdiabatic::~SurfCollideAdiabatic()
{
  if (copy) return;

  delete random;
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   ip = particle with current x = collision pt, current v = incident v
   isurf = index of surface element
   norm = surface normal unit vector
   isr = index of reaction model if >= 0, -1 for no chemistry
   ip = reset to NULL if destroyed by chemsitry
   return jp = new particle if created by chemistry
   return reaction = index of reaction (1 to N) that took place, 0 = no reaction
   resets particle(s) to post-collision outward velocity
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideAdiabatic::
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

  // note that adiabatic condition (i.e. not energy transfer of flow to surf)
  // does only apply to particle collisions. Chemistry (e.g. particle
  // adsorptions) can lead to energy transfer in both directions.
  if (isr >= 0) {
    if (modify->n_surf_react) memcpy(&iorig,ip,sizeof(Particle::OnePart));
    reaction = surf->sr[isr]->react(ip,isurf,norm,jp,velreset);
    if (reaction) surf->nreact_one++;
  }

  // isotropic scattering conserving velocity magnitude (i.e. kinetic energy)
  //   of each particle
  // only if SurfReact did not already reset velocities
  // cannot trigger fixes that require temperature of particle here
  //   because temperature of wall is not known

  if (ip) {
    if (!velreset) scatter_isotropic(ip,norm);
  }
  if (jp) {
    if (!velreset) scatter_isotropic(jp,norm);
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
   particle collision with adiabatic surface
   p = particle with current x = collision pt, current v = incident v
   norm = surface normal unit vector
   resets particle(s) to post-collision outward velocity so that particle
   is scattered isotropically whilst conserving its velocity magnitude
   (i.e. no energy transfer to surf), this is an efficient way to model
   an adiabatic surface in DSMC as it does not require an iterative
   procedure to determine the surface temp, for more details and references
   see documentation
------------------------------------------------------------------------- */

void SurfCollideAdiabatic::scatter_isotropic(Particle::OnePart *p, double *norm)
{
  double *v = p->v;
  double dot = MathExtra::dot3(v,norm);

  // tangent1/2 = surface tangential unit vectors

  double tangent1[3], tangent2[3];
  tangent1[0] = v[0] - dot*norm[0];
  tangent1[1] = v[1] - dot*norm[1];
  tangent1[2] = v[2] - dot*norm[2];

  // if mag(tangent1) == 0 mean normal collision, in that case choose
  // a random tangential vector
  if (MathExtra::lensq3(tangent1) == 0.0) {
    tangent2[0] = random->uniform();
    tangent2[1] = random->uniform();
    tangent2[2] = random->uniform();
    MathExtra::cross3(norm,tangent2,tangent1);
  }

  // normalize tangent1
  MathExtra::norm3(tangent1);
  // compute tangent2 as norm x tangent1, so that tangent1 and tangent2 are
  // orthogonal
  MathExtra::cross3(norm,tangent1,tangent2);

  // isotropic scattering
  // vmag = magnitude of incidient particle velocity vector
  // vperp = velocity component perpendicular to surface along norm (cy)
  // vtan1/2 = 2 remaining velocity components tangential to surface

  double vmag = MathExtra::len3(v);

  double theta = MY_2PI*random->uniform();
  double f_phi = random->uniform();
  double sqrt_f_phi = sqrt(f_phi);

  double vperp = vmag * sqrt(1.0 - f_phi);
  double vtan1 = vmag * sqrt_f_phi * sin(theta);
  double vtan2 = vmag * sqrt_f_phi * cos(theta);

  // update particle velocity
  v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
  v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
  v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

  // p->erot and p->evib stay identical
}

/* ----------------------------------------------------------------------
   wrapper on scatter_isotropic() method to perform collision for a
   single particle
   pass in 0 coefficients to match command-line args for style adiabatic
   flags, coeffs can be NULL
   called by SurfReactAdsorb
------------------------------------------------------------------------- */

void SurfCollideAdiabatic::wrapper(Particle::OnePart *p, double *norm,
                                  int *flags, double *coeffs)
{
  scatter_isotropic(p,norm);
}
