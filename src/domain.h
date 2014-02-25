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

#ifndef SPARTA_DOMAIN_H
#define SPARTA_DOMAIN_H

#include "pointers.h"
#include "particle.h"

namespace SPARTA_NS {

class Domain : protected Pointers {
 public:
  int box_exist;                    // 0 = not yet created, 1 = exists
  int dimension;                    // 2,3
  int axisymmetric;                 // 1 for yes, 0 for no, only allowed in 2d

  double boxlo[3],boxhi[3];         // box global bounds
  double xprd,yprd,zprd;            // global box dimensions
  double prd[3];                    // array form of dimensions

  int bflag[6];                     // boundary flags
  double norm[6][3];                // boundary normals

  int nregion;                      // # of defined Regions
  int maxregion;                    // max # regions can hold
  class Region **regions;           // list of defined Regions

  Domain(class SPARTA *);
  ~Domain();
  void init();
  void set_initial_box();
  void set_global_box();
  void set_boundary(int, char **);
  void boundary_modify(int, char **);
  int collide(Particle::OnePart *, int, int, double *);
  void uncollide(int, double *);
  void add_region(int, char **);
  void delete_region(int, char **);
  int find_region(char *);
  void print_box(const char *);

 private:
  int surf_collide[6];              // index of SurfCollide model
                                    // for each bflag = SURFACE boundary
};

}

#endif

/* ERROR/WARNING messages:

*/
