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

#ifndef DSMC_COMPUTE_H
#define DSMC_COMPUTE_H

#include "pointers.h"

namespace DSMC_NS {

class Compute : protected Pointers {
 public:
  char *id,*style;

  double scalar;            // computed global scalar
  double *vector;           // computed global vector
  double **array;           // computed global array
  double *vector_molecule;  // computed per-molecule vector
  double **array_molecule;  // computed per-molecule array
  double *vector_grid;      // computed per-grid vector
  double **array_grid;      // computed per-grid array
  double *vector_surf;      // computed per-surf vector
  double **array_surf;      // computed per-surf array

  int scalar_flag;          // 0/1 if compute_scalar() function exists
  int vector_flag;          // 0/1 if compute_vector() function exists
  int array_flag;           // 0/1 if compute_array() function exists
  int size_vector;          // length of global vector
  int size_array_rows;      // rows in global array
  int size_array_cols;      // columns in global array

  int per_molecule_flag;      // 0/1 if compute_per_molecule() function exists
  int size_per_molecule_cols; // 0 = vector, N = columns in per-molecule array

  int per_grid_flag;          // 0/1 if compute_per_grid() function exists
  int size_per_grid_cols;     // 0 = vector, N = columns in per-grid array

  int per_surf_flag;          // 0/1 if compute_per_surf() function exists
  int size_per_surf_cols;     // 0 = vector, N = columns in per-surf array

  int invoked_flag;       // non-zero if invoked or accessed this step, 0 if not
  bigint invoked_scalar;  // last timestep on which compute_scalar() was invoked
  bigint invoked_vector;       // ditto for compute_vector()
  bigint invoked_array;        // ditto for compute_array()
  bigint invoked_per_molecule; // ditto for compute_per_molecule()
  bigint invoked_per_grid;     // ditto for compute_per_grid()
  bigint invoked_per_surf;     // ditto for compute_per_surf()

  Compute(class DSMC *, int, char **);
  virtual ~Compute();
  virtual void init() {}

  virtual double compute_scalar() {return 0.0;}
  virtual void compute_vector() {}
  virtual void compute_array() {}
  virtual void compute_per_molecule() {}
  virtual void compute_per_grid() {}
  virtual void compute_per_surf() {}

  virtual double memory_usage() {return 0.0;}
};

}

#endif
