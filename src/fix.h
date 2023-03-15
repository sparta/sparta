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

#ifndef SPARTA_FIX_H
#define SPARTA_FIX_H

#include "pointers.h"
#include "particle.h"

namespace SPARTA_NS {

class Fix : protected Pointers {
 public:
  char *id,*style;

  int nevery;                    // how often to call an end_of_step fix
  int time_depend;               // 1 if requires continuous timestepping
  int gridmigrate;               // 0/1 if per grid cell info must migrate
  int flag_update_custom;         // 0/1 if has update_custom() method
  int flag_gas_react;            // 0/1 if has gas_react() method
  int flag_surf_react;           // 0/1 if has surf_react() method

  int scalar_flag;               // 0/1 if compute_scalar() function exists
  int vector_flag;               // 0/1 if compute_vector() function exists
  int array_flag;                // 0/1 if compute_array() function exists
  int size_vector;               // length of global vector
  int size_array_rows;           // rows in global array
  int size_array_cols;           // columns in global array
  int global_freq;               // frequency s/v data is available at

  int per_particle_flag;         // 0/1 if per-particle data is stored
  int size_per_particle_cols;    // 0 = vector, N = cols in per-particle array
  int per_particle_freq;         // frequency per-particle data is available at

  int per_grid_flag;             // 0/1 if per-grid data is stored
  int size_per_grid_cols;        // 0 = vector, N = cols in per-grid array
  int per_grid_freq;             // frequency per-grid data is available at

  int per_surf_flag;             // 0/1 if per-surf data is stored
  int size_per_surf_cols;        // 0 = vector, N = cols in per-surf array
  int per_surf_freq;             // frequency per-surf data is available at

  double *vector_particle;       // computed per-particle vector
  double **array_particle;       // computed per-particle array
  double *vector_grid;           // computed per-grid vector
  double **array_grid;           // computed per-grid array
  double *vector_surf;           // computed per-surf vector
  double **array_surf;           // computed per-surf array

  int per_particle_field;        // 0/1 if produces per-particle external field
  int per_grid_field;            // 0/1 if produces per-grid external field
  int field_active[3];           // 0/1 for active x,y,z components of ext field

  int START_OF_STEP,END_OF_STEP;    // mask settings

  int kokkos_flag;              // 0/1 if Kokkos fix
  int copy,copymode;            // 1 if copy of class (prevents deallocation of
                                //  base class when child copy is destroyed)
  ExecutionSpace execution_space;
  unsigned int datamask_read,datamask_modify;

  Fix(class SPARTA *, int, char **);
  Fix(class SPARTA *sparta) : Pointers(sparta) {} // needed for Kokkos
  virtual ~Fix();

  virtual int setmask() = 0;

  virtual void init() {}
  virtual void setup() {}

  virtual void start_of_step() {}
  virtual void end_of_step() {}
  virtual void update_custom(int, double, double, double, double *) {}
  virtual void gas_react(int) {}
  virtual void surf_react(Particle::OnePart *, int &, int &) {}
  virtual void compute_field() {}

  virtual int pack_grid_one(int, char *, int) {return 0;}
  virtual int unpack_grid_one(int, char *) {return 0;}
  virtual void copy_grid_one(int, int) {}
  virtual void add_grid_one() {}
  virtual void reset_grid_count(int) {}
  virtual void grid_changed() {}

  virtual double compute_scalar() {return 0.0;}
  virtual double compute_vector(int) {return 0.0;}
  virtual double compute_array(int,int) {return 0.0;}

  virtual double memory_usage() {return 0.0;}
};

}

#endif

/* ERROR/WARNING messages:

E: Fix ID must be alphanumeric or underscore characters

Self-explanatory.

*/
