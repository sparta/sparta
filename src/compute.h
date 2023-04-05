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

#ifndef SPARTA_COMPUTE_H
#define SPARTA_COMPUTE_H

#include "pointers.h"
#include "particle.h"

namespace SPARTA_NS {

class Compute : protected Pointers {
 public:
  char *id,*style;

  double scalar;            // computed global scalar
  double *vector;           // computed global vector
  double **array;           // computed global array
  double *vector_particle;  // computed per-particle vector
  double **array_particle;  // computed per-particle array
  double *vector_grid;      // computed per-grid vector
  double **array_grid;      // computed per-grid array

  // vec/array surf are length = # of explicit surf elements owned
  // vec/array tally are length = # of surf elements tallied

  double *vector_surf;        // computed per-surf vector
  double **array_surf;        // computed per-surf array
  double *vector_surf_tally;  // computed per-surf tally vector
  double **array_surf_tally;  // computed per-surf tally array

  int scalar_flag;          // 0/1 if compute_scalar() function exists
  int vector_flag;          // 0/1 if compute_vector() function exists
  int array_flag;           // 0/1 if compute_array() function exists
  int size_vector;          // length of global vector
  int size_array_rows;      // rows in global array
  int size_array_cols;      // columns in global array

  int per_particle_flag;      // 0/1 if compute_per_particle() function exists
  int size_per_particle_cols; // 0 = vector, N = columns in per-particle array

  int per_grid_flag;            // 0/1 if compute_per_grid() function exists
  int size_per_grid_cols;       // 0 = vector, N = columns in per-grid array
  int post_process_grid_flag;   // 1 if requires post_process_grid() for output
  int post_process_isurf_grid_flag; // 1 if requires post_process_tally() for out

  int per_surf_flag;          // 0/1 if compute_per_surf() function exists
  int size_per_surf_cols;     // 0 = vector, N = columns in per-surf array

  int surf_tally_flag;        // 1 if compute tallies surface bounce info
  int boundary_tally_flag;    // 1 if compute tallies boundary bounce info

  int timeflag;       // 1 if Compute stores list of timesteps it's called on
  int ntime;          // # of entries in time list
  int maxtime;        // max # of entries time list can hold
  bigint *tlist;      // list of timesteps the Compute is called on

  int invoked_flag;       // non-zero if invoked or accessed this step, 0 if not
  bigint invoked_scalar;  // last timestep on which compute_scalar() was invoked
  bigint invoked_vector;       // ditto for compute_vector()
  bigint invoked_array;        // ditto for compute_array()
  bigint invoked_per_particle; // ditto for compute_per_particle()
  bigint invoked_per_grid;     // ditto for compute_per_grid()
  bigint invoked_per_surf;     // ditto for compute_per_surf()

  Compute(class SPARTA *, int, char **);
  Compute(class SPARTA* sparta) : Pointers(sparta) {}
  virtual ~Compute();
  virtual void init() {}

  virtual double compute_scalar() {return 0.0;}
  virtual void compute_vector() {}
  virtual void compute_array() {}
  virtual void compute_per_particle() {}
  virtual void compute_per_grid() {}
  virtual void compute_per_surf() {}
  virtual void clear() {}
  virtual void surf_tally(int, int, int, Particle::OnePart *,
                          Particle::OnePart *, Particle::OnePart *) {}
  virtual void boundary_tally(int, int, int, Particle::OnePart *,
                              Particle::OnePart *, Particle::OnePart *) {}

  virtual void post_process_grid(int, int, double **, int *, double *, int) {}
  // NOTE: get rid of this method at some point
  virtual void post_process_grid_old(void *, void *, int, int, double *, int) {}
  virtual void post_process_isurf_grid() {}

  virtual int query_tally_grid(int, double **&, int *&) {return 0;}
  virtual int tallyinfo(surfint *&) {return 0;}
  virtual void post_process_surf() {}

  virtual void reallocate() {}
  virtual bigint memory_usage();

  // methods in compute.cpp

  void addstep(bigint);
  int matchstep(bigint);
  void clearstep();

  // Kokkos methods

  int kokkos_flag;          // 1 if Kokkos-enabled
  int copy,copymode;        // 1 if copy of class (prevents deallocation of
                            //  base class when child copy is destroyed)
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute ID must be alphanumeric or underscore characters

Self-explanatory.

*/
