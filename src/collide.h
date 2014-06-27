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
  int rotstyle;       // none/smooth rotational modes
  int vibstyle;       // none/discrete/smooth vibrational modes

  int ncollide_one,nattempt_one,nreact_one;
  bigint ncollide_running,nattempt_running,nreact_running;
 
  Collide(class SPARTA *, int, char **);
  virtual ~Collide();
  virtual void init();
  void modify_params(int, char **);
  void reset_vremax();
  void collisions();

  virtual double vremax_init(int, int) = 0;
  virtual double attempt_collision(int, int, double) = 0;
  virtual double attempt_collision(int, int, int, double) = 0;
  virtual int test_collision(int, int, int, 
			     Particle::OnePart *, Particle::OnePart *) = 0;
  virtual void setup_collision(Particle::OnePart *, Particle::OnePart *) = 0;
  virtual Particle::OnePart *perform_collision(Particle::OnePart *, 
					       Particle::OnePart *) = 0;

  virtual double extract(int, const char *) {return 0.0;}

  int pack_grid_one(int, char *, int);
  int unpack_grid_one(int, char *);
  void compress_grid();

 protected:
  int npmax;          // max # of particles in plist
  int *plist;         // list of particles in a single group

  int nglocal;        // current size of per-cell arrays
  int nglocalmax;     // max allocated size of per-cell arrays

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

  int vre_first;      // 1 for first run after collision style is defined
  int vre_start;      // 1 if reset vre params at start of each run
  int vre_every;      // reset vre params every this many steps
  bigint vre_next;    // next timestep to reset vre params on
  int remainflag;     // 1 if remain defined, else use random fraction

  double ***vremax;   // max relative velocity, per cell, per group pair
  double ***remain;   // collision number remainder, per cell, per group pair
  double **vremax_initial;   // initial vremax value, per group pair

  inline void addgroup(int igroup, int n)
  {
    if (ngroup[igroup] == maxgroup[igroup]) {
      maxgroup[igroup] += DELTAPART;
      memory->grow(glist[igroup],maxgroup[igroup],"collide:grouplist");
    }
    glist[igroup][ngroup[igroup]++] = n;
  }

  void collisions_one();
  void collisions_group();
  void grow_percell(int);
};

}

#endif

/* ERROR/WARNING messages:

E: Collision mixture does not exist

Self-explantory.

E: Collision mixture does not contain all species

The specified mixture must contain all species in the simulation so
that they can be assigned to collision groups.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/
