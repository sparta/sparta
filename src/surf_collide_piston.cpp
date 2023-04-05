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
#include "surf_collide_piston.h"
#include "surf.h"
#include "surf_react.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "input.h"
#include "math_extra.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

SurfCollidePiston::SurfCollidePiston(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal surf_collide piston command");

  vwall = input->numeric(FLERR,arg[2]);
  if (vwall <= 0.0) error->all(FLERR,"Surf_collide piston velocity <= 0.0");
}

/* ---------------------------------------------------------------------- */

void SurfCollidePiston::init()
{
  SurfCollide::init();

  dt = update->dt;

  // check that this model only assigned to surfs with axis-aligned normals
  // index = position in surf->sc list that this SurfCollide instance is

  int index;
  for (index = 0; index < surf->nsc; index++)
    if (this == surf->sc[index]) break;

  int flag = 0;

  if (domain->dimension == 2) {
    Surf::Line *lines = surf->lines;
    int nsurf = surf->nsurf;
    for (int i = 0; i < nsurf; i++)
      if (lines[i].isc == index) {
        if (lines[i].norm[0] != 0.0 && lines[i].norm[1] != 0.0) flag++;
      }
  }

  if (domain->dimension == 3) {
    Surf::Tri *tris = surf->tris;
    int nsurf = surf->nsurf;
    for (int i = 0; i < nsurf; i++)
      if (tris[i].isc == index) {
        if (tris[i].norm[0] != 0.0 && tris[i].norm[1] != 0.0) flag++;
        if (tris[i].norm[1] != 0.0 && tris[i].norm[2] != 0.0) flag++;
        if (tris[i].norm[2] != 0.0 && tris[i].norm[0] != 0.0) flag++;
      }
  }

  if (flag) error->all(FLERR,"Surf_collide piston assigned to "
                       "surface with non axis-aligned normal");
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   ip = particle with current x = collision pt, current v = incident v
   dtremain = portion of timestep remaining
   isurf = index of surface element
   norm = surface normal unit vector
   isr = index of reaction model if >= 0, -1 for no chemistry
   ip = reset to NULL if destroyed by chemistry
   return jp = new particle if created by chemistry
   return reaction = index of reaction (1 to N) that took place, 0 = no reaction
   resets particle(s) to post-collision outward velocity
   update dtremain
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollidePiston::
collide(Particle::OnePart *&ip, double &dtremain,
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

  // norm will be in single coordinate direction
  // dir = 0,1,2 for wall (or surface) with norm parallel to x,y,z
  // which = 0/1 for wall (or surface) with +/- normal (lo/hi wall)

  int dim,which;

  if (norm[0] != 0.0) {
    dim = 0;
    if (norm[0] < 0.0) which = 1;
    else which = 0;
  } else if (norm[1] != 0.0) {
    dim = 1;
    if (norm[1] < 0.0) which = 1;
    else which = 0;
  } else {
    dim = 2;
    if (norm[2] < 0.0) which = 1;
    else which = 0;
  }

  // xwall = initial position of wall (collision pt)
  // xorig = initial coordinate component
  // vorig = initial velocity component
  // vwall = user-specified wall velocity (always >= 0)

  double *x = ip->x;
  double *v = ip->v;
  double xwall = x[dim];
  double xorig = xwall - v[dim]*(dt - dtremain);
  double vorig = v[dim];

  // piston reflection: see eqs 12.30 and 12.31 in Bird 1994, p 288
  // uprime = post-collision velocity component
  // xprime = post-collision coordinate component
  // delete particle and return if xprime is not inside box
  // formula for dtremain works for both which = 0/1
  //   since numerator and denominator are always same sign

  double uprime,xprime;

  if (which == 0) {
    uprime = -2.0*vwall - vorig;
    xprime = 2.0*xwall - xorig + uprime*dt;
    if (xprime <= xwall) {
      ip = NULL;
      return NULL;
    }
  } else {
    uprime = 2.0*vwall - vorig;
    xprime = 2.0*xwall - xorig + uprime*dt;
    if (xprime >= xwall) {
      ip = NULL;
      return NULL;
    }
  }

  dtremain = (xprime - xwall) / uprime;
  v[dim] = uprime;

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
