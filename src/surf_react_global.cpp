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
#include "surf_react_global.h"
#include "input.h"
#include "update.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_extra.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

SurfReactGlobal::SurfReactGlobal(SPARTA *sparta, int narg, char **arg) :
  SurfReact(sparta, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal surf_react global command");

  prob_destroy = input->numeric(FLERR,arg[2]);
  prob_create = input->numeric(FLERR,arg[3]);

  if (prob_destroy + prob_create > 1.0)
    error->all(FLERR,"Illegal surf_react global command");

  // setup the reaction tallies

  nsingle = ntotal = 0;

  nlist = 2;
  tally_single = new int[nlist];
  tally_total = new int[nlist];
  tally_single_all = new int[nlist];
  tally_total_all = new int[nlist];

  size_vector = 2 + 2*nlist;

  // initialize RNG

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfReactGlobal::~SurfReactGlobal()
{
  if (copy) return;

  delete random;
}

/* ----------------------------------------------------------------------
   select surface reaction to perform for particle with ptr IP on surface
   return which reaction 1 (destroy), 2 (create), 0 = no reaction
   if create, add particle and return ptr JP
------------------------------------------------------------------------- */

int SurfReactGlobal::react(Particle::OnePart *&ip, int, double *,
                           Particle::OnePart *&jp, int &)
{
  double r = random->uniform();

  // perform destroy reaction

  if (r < prob_destroy) {
    nsingle++;
    tally_single[0]++;
    ip = NULL;
    return 1;
  }

  // perform create reaction
  // clone 1st particle to create 2nd particle
  // if add_particle performs a realloc:
  //   make copy of x,v with new species
  //   rot/vib energies will be reset by SurfCollide
  //   repoint ip to new particles data struct if reallocated

  if (r < prob_destroy+prob_create) {
    nsingle++;
    tally_single[1]++;
    double x[3],v[3];
    int id = MAXSMALLINT*random->uniform();
    memcpy(x,ip->x,3*sizeof(double));
    memcpy(v,ip->v,3*sizeof(double));
    Particle::OnePart *particles = particle->particles;
    int reallocflag =
      particle->add_particle(id,ip->ispecies,ip->icell,x,v,0.0,0.0);
    if (reallocflag) ip = particle->particles + (ip - particles);
    jp = &particle->particles[particle->nlocal-1];
    return 2;
  }

  // no reaction

  return 0;
}

/* ---------------------------------------------------------------------- */

char *SurfReactGlobal::reactionID(int m)
{
  if (m == 0) return (char *) "delete";
  return (char *) "create";
}
