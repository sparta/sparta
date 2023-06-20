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

#ifndef SPARTA_SURF_COLLIDE_H
#define SPARTA_SURF_COLLIDE_H

#include "pointers.h"
#include "particle.h"

namespace SPARTA_NS {

class SurfCollide : protected Pointers {
 public:
  char *id;
  char *style;

  int allowreact;           // 1 if allows for surface reactions
  int dynamicflag;          // 1 if any param is dynamically updated
  int transparent;          // 1 if transparent collision model
  int vector_flag;          // 0/1 if compute_vector() function exists
  int size_vector;          // length of global vector

  SurfCollide(class SPARTA *, int, char **);
  SurfCollide(class SPARTA *sparta) : Pointers(sparta) {}
  virtual ~SurfCollide();
  virtual void init();
  virtual Particle::OnePart *collide(Particle::OnePart *&, double &,
                                     int, double *, int, int &) = 0;
  virtual void wrapper(Particle::OnePart *, double *, int *, double *) {}
  virtual void flags_and_coeffs(int *, double *) {}

  void dynamic();
  void tally_reset();
  void tally_update();
  double compute_vector(int i);

  int copy,copymode;

 protected:

  // tallies for collisions
  // nsingle = all collisions in one step
  // ntotal = cumulative nsingle across all steps
  // one,all used in compute_vector()

  int nsingle,ntotal;
  double one[2],all[2];

  // variables used by all SC classes which define Tsurf
  
  int tmode;                 // possible modes = NUMERIC,VAREQUAL,VARSURF,CUSTOM
  double tsurf;
  char *tstr;                // temperature variable name (NULL if constant)
  int tindex_var;                  // index of equal-style variable
  int tindex_custom;                  // index of equal-style variable
  double *tvector;           // custom per-surf temperature vector
  int tfreq;                 // frequency to update variables
  int persurf_temperature;
  double *t_persurf;
  
  // functions used by all SC classes which define Tsurf
  
  void parse_tsurf(char *);
  void check_tsurf();
};

}

#endif

/* ERROR/WARNING messages:

E: Surf_collide ID must be alphanumeric or underscore characters

Self-explanatory.

*/
