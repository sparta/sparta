/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_UPDATE_H
#define SPARTA_UPDATE_H

#include "pointers.h"

namespace SPARTA_NS {

class Update : protected Pointers {
 public:
  bigint ntimestep;               // current timestep
  int nsteps;                     // # of steps to run
  bigint firststep,laststep;      // 1st & last step of this run
  int runflag;                    // 0 for unset, 1 for run
  double dt;                      // timestep size

  char *unit_style;      // style of units used throughout simulation
  double boltz;          // Boltzmann constant (eng/degree K)
  double mvv2e;          // conversion of mv^2 to energy

  double fnum;           // ratio of real particles to simulation particles
  double nrho;           // number density of background gas
  double vstream[3];     // streaming velocity of background gas
  double temp_thermal;   // thermal temperature of background gas

  int nmigrate;          // # of particles to migrate to new procs
  int *mlist;            // indices of particles to migrate

                         // current step counters
  int ntouch_one;        // particle-cell touches
  int ncomm_one;         // particles migrating to new procs
  int nboundary_one;     // particles colliding with global boundary
  int nexit_one;         // particles exiting outflow boundary
  int nscheck_one;       // surface elements checked for collisions
  int nscollide_one;     // particle/surface collisions

  bigint nmove_running;      // running count of total particle moves
  bigint ntouch_running;     // running count of current step counters
  bigint ncomm_running;
  bigint nboundary_running;
  bigint nexit_running;
  bigint nscheck_running;
  bigint nscollide_running;

  class RanMars *ranmaster;   // master random number generator

  Update(class SPARTA *);
  ~Update();
  void init();
  void set_units(const char *);
  void setup();
  void run(int);
  void global(int, char **);

 private:
  int me,nprocs;
  int ncurrent;              // local # of particles before insertion
  int maxmigrate;            // max # of particles in mlist
  int faceflip[6];

  class RanPark *random;     // RNG for particle timestep moves

  int bounce_tally;                  // 1 if any bounces tallied on this step
  int surf_tally_flag;               // 1 if tally surf bounces on this step
  int boundary_tally_flag;           // 1 if tally boundary bounces on this step

  int nslist_compute;                // # of surf bounce computes to check
  int nblist_compute;                // # of boundary bounce computes to check
  class Compute **slist_compute;     // list of surf bounce Computes
  class Compute **blist_compute;     // list of boundary bounce Computes
  
  int bounce_setup();
  void bounce_set(bigint);

  typedef void (Update::*FnPtr)();
  FnPtr move;                // ptr to move method

  void move3d_surface();     // variants of move method
  void move3d();
  void move2d_surface();
  void move2d();
};

}

#endif
