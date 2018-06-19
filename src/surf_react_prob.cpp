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
#include "surf_react_prob.h"
#include "input.h"
#include "update.h"
#include "comm.h"
#include "random_mars.h"
#include "random_park.h"
#include "math_extra.h"
#include "error.h"

using namespace SPARTA_NS;

enum{DISSOCIATION,EXCHANGE,RECOMBINATION};        // other surf react files

/* ---------------------------------------------------------------------- */

SurfReactProb::SurfReactProb(SPARTA *sparta, int narg, char **arg) :
  SurfReact(sparta, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal surf_react prob command");

  readfile(arg[2]);

  // initialize RNG

  random = new RanPark(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfReactProb::~SurfReactProb()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

void SurfReactProb::init()
{
  SurfReact::init();
  init_reactions();
}

/* ---------------------------------------------------------------------- */

int SurfReactProb::react(Particle::OnePart *&ip, double *, 
                         Particle::OnePart *&jp)
{
  int n = reactions[ip->ispecies].n;

  if (n == 0) return 0;
  int *list = reactions[ip->ispecies].list;

  // probablity to compare to reaction probability

  double react_prob = 0.0;
  double random_prob = random->uniform(); 

  // loop over possible reactions for this species
  // if dissociation performs a realloc:
  //   make copy of x,v with new species
  //   rot/vib energies will be reset by SurfCollide
  //   repoint ip to new particles data struct if reallocated

  OneReaction *r;

  for (int i = 0; i < n; i++) {
    r = &rlist[list[i]];
    react_prob += r->coeff[0];

    if (react_prob > random_prob) {
      nsingle++;
      switch (r->type) {
      case DISSOCIATION:
        {
          double x[3],v[3];
          ip->ispecies = r->products[0];
          int id = MAXSMALLINT*random->uniform();
          memcpy(x,ip->x,3*sizeof(double));
          memcpy(v,ip->v,3*sizeof(double));  
          Particle::OnePart *particles = particle->particles;
          int reallocflag = 
            particle->add_particle(id,r->products[1],ip->icell,x,v,0.0,0.0);
          if (reallocflag) ip = particle->particles + (ip - particles);
          jp = &particle->particles[particle->nlocal-1];
          return 1;
        }
      case EXCHANGE:
        {
          ip->ispecies = r->products[0];
          return 1;
        }
      case RECOMBINATION:
        {
          ip = NULL;
          return 1;
        }
      }
    }
  }

  // no reaction

  return 0;
}
