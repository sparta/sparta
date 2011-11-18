/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#ifndef DSMC_COLLIDE_H
#define DSMC_COLLIDE_H

#include "pointers.h"
#include "particle.h"

namespace DSMC_NS {

class Collide : protected Pointers {
 public:
  bigint ncollattempt;
  bigint ncollision;
  char *style;
 
  Collide(class DSMC *, int, char **);
  virtual ~Collide();
  virtual void init();
  void collisions();

  virtual int attempt_collision(int, int, int, double) = 0;
  virtual int test_collision(int, int, int, 
			     Particle::OnePart *, Particle::OnePart *) = 0;
  virtual void setup_collision(Particle::OnePart *, Particle::OnePart *) = 0;
  virtual Particle::OnePart *perform_collision(Particle::OnePart *, 
					       Particle::OnePart *) = 0;

 protected:
  int nspecies;       // # of species defined
  int oldspecies;     // # of species defined on previous run
  int *nsp;           // # of particles in one cell of each species
  int *maxsp;         // max # of splist indices allocated per species
  int **splist;       // indices of particles in one cell of each species

  int nsspair;        // # of species pairs to do collisions for
  int **sscoll;       // Nx3 list of species pairs to do collisions for
                      // 0 = sp I, 1 = sp J, 2 = # of collisions to attempt

  class RanPark *random;
};

}

#endif
