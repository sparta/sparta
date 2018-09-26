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
#include <iostream>

namespace SPARTA_NS {

enum{COPYPARTICLELIST,FIXEDMEMORY};

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

  int reorder_period;        // # of timesteps between particle reordering
  int nParticlesReorderSet;  // # of particles in reorder set in fixed memory reordering scheme
                             //  (on GPU, this will be number of threads for each reordering pass)
  int reorder_scheme;        // (copying full particle list or using a small fixed memory)

  int copymode;          // 1 if copy of class (prevents deallocation of
                         //  base class when child copy is destroyed)

  class RanMars *ranmaster;   // master random number generator

  double rcblo[3],rcbhi[3];    // debug info from RCB for dump image

  Update(class SPARTA *);
  ~Update();
  void set_units(const char *);
  virtual void init();
  virtual void setup();
  virtual void run(int);
  void global(int, char **);
  void reset_timestep(int, char **);

  int split3d(int, double *);
  int split2d(int, double *);

 protected:
  int me,nprocs;
  int maxmigrate;            // max # of particles in mlist
  class RanPark *random;     // RNG for particle timestep moves

  int collide_react;         // 1 if any SurfCollide or React classes defined
  int nsc,nsr;               // copy of Collide/React data in Surf class
  class SurfCollide **sc;
  class SurfReact **sr;

  int bounce_tally;               // 1 if any bounces are ever tallied
  int nslist_compute;             // # of computes that tally surf bounces
  int nblist_compute;             // # of computes that tally boundary bounces
  class Compute **slist_compute;  // list of all surf bounce Computes
  class Compute **blist_compute;  // list of all boundary bounce Computes

  int nsurf_tally;         // # of Cmp tallying surf bounce info this step
  int nboundary_tally;     // # of Cmp tallying boundary bounce info this step
  class Compute **slist_active;   // list of active surf Computes this step
  class Compute **blist_active;   // list of active boundary Computes this step

  int surf_pre_tally;       // 1 to log particle stats before surf collide
  int boundary_pre_tally;   // 1 to log particle stats before boundary collide

  int collide_react_setup();
  void collide_react_update();

  int bounce_setup();
  virtual void bounce_set(bigint);
  void reset_timestep(bigint);

  //int axi_vertical_line(double, double *, double *, double, double, double,
  //                     double &);

  // remap x and v components into axisymmetric plane
  // input x at end of linear move (x = xold + dt*v)
  // change x[1] = sqrt(x[1]^2 + x[2]^2), x[2] = 0.0
  // change vy,vz by rotation into axisymmetric plane

  inline void axi_remap(double *x, double *v) {
    double ynew = x[1];
    double znew = x[2];
    x[1] = sqrt(ynew*ynew + znew*znew);
    x[2] = 0.0;
    double rn = ynew / x[1];
    double wn = znew / x[1];
    double vy = v[1];
    double vz = v[2];
    v[1] = vy*rn + vz*wn;
    v[2] = -vy*wn + vz*rn;
  };

  typedef void (Update::*FnPtr)();
  FnPtr moveptr;             // ptr to move method
  template < int, int > void move();

  int perturbflag;
  typedef void (Update::*FnPtr2)(double, double *, double *);
  FnPtr2 moveperturb;        // ptr to moveperturb method

  // variants of moveperturb method
  // adjust end-of-move x,v due to perturbation on straight-line advection

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

Self-explanatory.

E: Gravity in y not allowed for axi-symmetric model

Self-explanatory.

E: Particle %d on proc %d hit inside of surf %d on step %ld

This error should not happen if particles start outside of physical
objects.  Please report the issue to the SPARTA developers.

E: Sending particle to self

This error should not occur.  Please report the issue to the SPARTA
developers.

E: Cannot set global surfmax when surfaces already exist

This setting must be made before any surfac elements are
read via the read_surf command.

E: Timestep must be >= 0

Reset_timestep cannot be used to set a negative timestep.

E: Too big a timestep

Reset_timestep timestep value must fit in a SPARTA big integer, as
specified by the -DSPARTA_SMALL, -DSPARTA_BIG, or -DSPARTA_BIGBIG
options in the low-level Makefile used to build SPARTA.  See
Section 2.2 of the manual for details.

E: Cannot reset timestep with a time-dependent fix defined

The timestep cannot be reset when a fix that keeps track of elapsed
time is in place.

*/
