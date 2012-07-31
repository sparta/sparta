/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   www.sandia.gov/sparta.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_COLLIDE_H
#define SPARTA_COLLIDE_H

#include "pointers.h"
#include "memory.h"
#include "particle.h"

namespace SPARTA_NS {

#define DELTAPART 128

class Collide : protected Pointers {
 public:
  char *style;

  int ncollide_one,nattempt_one;
  bigint ncollide_running,nattempt_running;
 
  Collide(class SPARTA *, int, char **);
  virtual ~Collide();
  virtual void init();
  void collisions();

  virtual double attempt_collision(int, int, int, double) = 0;
  virtual int test_collision(int, int, int, 
			     Particle::OnePart *, Particle::OnePart *) = 0;
  virtual void setup_collision(Particle::OnePart *, Particle::OnePart *) = 0;
  virtual Particle::OnePart *perform_collision(Particle::OnePart *, 
					       Particle::OnePart *) = 0;

 protected:
  int ngroups;        // # of groups
  int *ngroup;        // # of particles in one cell of each group
  int *maxgroup;      // max # of glist indices allocated per group
  int **glist;        // indices of particles in one cell of each group

  int npair;          // # of group pairs to do collisions for
  int **gpair;        // Nx3 list of species pairs to do collisions for
                      // 0 = igroup, 1 = jgroup, 2 = # of attempt collisions

  char *mixID;               // ID of mixture to use for groups
  class Mixture *mixture;    // ptr to mixture
  class RanPark *random;     // RNG for collision generation

  inline void addgroup(int igroup, int n)
  {
    if (ngroup[igroup] == maxgroup[igroup]) {
      maxgroup[igroup] += DELTAPART;
      memory->grow(glist[igroup],maxgroup[igroup],"collide:grouplist");
    }
    glist[igroup][ngroup[igroup]++] = n;
  }
};

}

#endif
