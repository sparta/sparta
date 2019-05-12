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
#include "ctype.h"
#include "stdlib.h"
#include "string.h"
#include "surf_react_zuzax.h"
#include "input.h"
#include "update.h"
#include "comm.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"


enum{DISSOCIATION,EXCHANGE,RECOMBINATION};        // other surf react files
enum{SIMPLE};                                     // other surf react files

#define MAXREACTANT 1
#define MAXPRODUCT 2
#define MAXCOEFF 2

#define MAXLINE 1024
#define DELTALIST 16

namespace SPARTA_NS {


//=================================================================================================
SurfReactZuzax::SurfReactZuzax(SPARTA *sparta, int narg, char **arg) : 
    SurfReact(sparta, narg, arg)
{

}

//=================================================================================================

SurfReactZuzax::~SurfReactZuzax()
{
}

/* ---------------------------------------------------------------------- */

void SurfReactZuzax::init()
{
  nsingle = ntotal = 0;
}

/* ---------------------------------------------------------------------- */

void SurfReactZuzax::init_reactions() 
{
}
/* ---------------------------------------------------------------------- */

int SurfReactZuzax::react(Particle::OnePart *&ip, double *tmpp, Particle::OnePart *&jp)
{
  double r = random->uniform();

  // perform destroy reaction

  if (r < prob_destroy) {
    nsingle++;
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
    double x[3],v[3];
    int id = MAXSMALLINT*random->uniform();
    memcpy(x,ip->x,3*sizeof(double));
    memcpy(v,ip->v,3*sizeof(double));
    Particle::OnePart *particles = particle->particles;
    int reallocflag =
      particle->add_particle(id,ip->ispecies,ip->icell,x,v,0.0,0.0);
    if (reallocflag) ip = particle->particles + (ip - particles);
    jp = &particle->particles[particle->nlocal-1];
    return 1;
  }

  // no reaction

  return 0;
}
//=================================================================================================

