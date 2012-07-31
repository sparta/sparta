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

#ifndef SPARTA_FIX_H
#define SPARTA_FIX_H

#include "pointers.h"

namespace SPARTA_NS {

class Fix : protected Pointers {
 public:
  char *id,*style;

  int scalar_flag;               // 0/1 if compute_scalar() function exists
  int vector_flag;               // 0/1 if compute_vector() function exists
  int array_flag;                // 0/1 if compute_array() function exists
  int size_vector;               // length of global vector
  int size_array_rows;           // rows in global array
  int size_array_cols;           // columns in global array
  int global_freq;               // frequency s/v data is available at

  int per_molecule_flag;         // 0/1 if per-molecule data is stored
  int size_per_molecule_cols;    // 0 = vector, N = cols in per-molecule array
  int per_molecule_freq;         // frequency per-molecule data is available at

  int per_grid_flag;             // 0/1 if per-grid data is stored
  int size_per_grid_cols;        // 0 = vector, N = cols in per-grid array
  int per_grid_freq;             // frequency per-grid data is available at

  int per_surf_flag;             // 0/1 if per-surf data is stored
  int size_per_surf_cols;        // 0 = vector, N = cols in per-surf array
  int per_surf_freq;             // frequency per-surf data is available at

  double *vector_molecule;       // computed per-molecule vector
  double **array_molecule;       // computed per-molecule array
  double *vector_grid;           // computed per-grid vector
  double **array_grid;           // computed per-grid array
  double *vector_surf;           // computed per-surf vector
  double **array_surf;           // computed per-surf array

  int bounceflag;                // 1 if fix needs surface bounce info
  bigint bouncenext;             // next timestep fix needs bounce info

  int START_OF_STEP,END_OF_STEP;    // mask settings

  Fix(class SPARTA *, int, char **);
  virtual ~Fix();

  virtual int setmask() = 0;

  virtual void init() {}
  virtual void setup(int) {}

  virtual void start_of_step() {}
  virtual void end_of_step() {}

  virtual double compute_scalar() {return 0.0;}
  virtual double compute_vector(int) {return 0.0;}
  virtual double compute_array(int,int) {return 0.0;}

  virtual double memory_usage() {return 0.0;}
};

}

#endif
