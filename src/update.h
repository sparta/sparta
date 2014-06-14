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

#ifndef SPARTA_UPDATE_H
#define SPARTA_UPDATE_H

#include "math.h"
#include "pointers.h"

namespace SPARTA_NS {

class Update : protected Pointers {
 public:
  bigint ntimestep;               // current timestep
  int nsteps;                     // # of steps to run
  int runflag;                    // 0 for unset, 1 for run
  bigint firststep,laststep;      // 1st & last step of this run
  bigint beginstep,endstep;       // 1st and last step of multiple runs
  int first_update;               // 0 before initial update, 1 after
  double dt;                      // timestep size

  char *unit_style;      // style of units used throughout simulation
  double boltz;          // Boltzmann constant (eng/degree K)
  double mvv2e;          // conversion of mv^2 to energy

  double fnum;           // ratio of real particles to simulation particles
  double nrho;           // number density of background gas
  double vstream[3];     // streaming velocity of background gas
  double temp_thermal;   // thermal temperature of background gas
  double gravity[3];     // acceleration vector of gravity

  int nmigrate;          // # of particles to migrate to new procs
  int *mlist;            // indices of particles to migrate

                         // current step counters
  int niterate;          // iterations of move/comm
  int ntouch_one;        // particle-cell touches
  int ncomm_one;         // particles migrating to new procs
  int nboundary_one;     // particles colliding with global boundary
  int nexit_one;         // particles exiting outflow boundary
  int nscheck_one;       // surface elements checked for collisions
  int nscollide_one;     // particle/surface collisions

  bigint first_running_step; // timestep running counts start on
  int niterate_running;      // running count of move/comm interations
  bigint nmove_running;      // running count of total particle moves
  bigint ntouch_running;     // running count of current step counters
  bigint ncomm_running;
  bigint nboundary_running;
  bigint nexit_running;
  bigint nscheck_running;
  bigint nscollide_running;

  int nstuck;                // # of particles stuck on surfs and deleted

  class RanMars *ranmaster;   // master random number generator

  double rcblo[3],rcbhi[3];    // debug info from RCB for dump image

  Update(class SPARTA *);
  ~Update();
  void set_units(const char *);
  void init();
  void setup();
  void run(int);
  void global(int, char **);
  void reset_timestep(int, char **);

  int split3d(int, double *);
  int split2d(int, double *);

 private:
  int me,nprocs;
  int maxmigrate;            // max # of particles in mlist
  class RanPark *random;     // RNG for particle timestep moves

  int bounce_tally;               // 1 if any bounces are ever tallied
  int nslist_compute;             // # of computes that tally surf bounces
  int nblist_compute;             // # of computes that tally boundary bounces
  class Compute **slist_compute;  // list of all surf bounce Computes
  class Compute **blist_compute;  // list of all boundary bounce Computes

  int nsurf_tally;         // # of Cmp tallying surf bounce info this step
  int nboundary_tally;     // # of Cmp tallying boundary bounce info this step
  class Compute **slist_active;   // list of active surf Computes this step
  class Compute **blist_active;   // list of active boundary Computes this step
  
  int bounce_setup();
  void bounce_set(bigint);
  void reset_timestep(bigint);

  typedef void (Update::*FnPtr)();
  FnPtr moveptr;             // ptr to move method
  template < int, int > void move();

  int perturbflag;
  typedef void (Update::*FnPtr2)(double, double *, double *);
  FnPtr2 moveperturb;        // ptr to moveperturb method

  // variants of moveperturb method

  inline void axisymmetric(double dt, double *x, double *v) {
    double dz = dt*v[2];
    double rold = x[1];
    x[1] = sqrt(x[1]*x[1] + dz*dz);
    double rn = rold / x[1];
    double wn = dz / x[1];
    double vold = v[1];
    double wold = v[2];
    v[1] = vold*rn + wold*wn;
    v[2] = -vold*wn + wold*rn;
  };

  inline void gravity2d(double dt, double *x, double *v) {
    double dtsq = 0.5*dt*dt;
    x[0] += dtsq*gravity[0];
    x[1] += dtsq*gravity[1];
    v[0] += dt*gravity[0];
    v[1] += dt*gravity[1];
  };

  inline void gravity3d(double dt, double *x, double *v) {
    double dtsq = 0.5*dt*dt;
    x[0] += dtsq*gravity[0];
    x[1] += dtsq*gravity[1];
    x[2] += dtsq*gravity[2];
    v[0] += dt*gravity[0];
    v[1] += dt*gravity[1];
    v[2] += dt*gravity[2];
  };
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Gravity in z not allowed for 2d

UNDOCUMENTED

E: Gravity in y not allowed for axi-symmetry

UNDOCUMENTED

E: Particle %d on proc %d hit inside of surf %d on step %ld

UNDOCUMENTED

E: Sending particle to self

UNDOCUMENTED

E: Cannot set global surfmax when surfaces already exist

UNDOCUMENTED

E: Timestep must be >= 0

UNDOCUMENTED

E: Too big a timestep

UNDOCUMENTED

E: Cannot reset timestep with a time-dependent fix defined

UNDOCUMENTED

*/
