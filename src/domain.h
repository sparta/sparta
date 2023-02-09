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
  int boundary_collision_check;  // flag for whether init() check is required
                                 // for assign of collision models to boundaries

  double boxlo[3],boxhi[3];         // box global bounds
  double xprd,yprd,zprd;            // global box dimensions
  double prd[3];                    // array form of dimensions

  int bflag[6];                     // boundary flags
  double norm[6][3];                // boundary normals

  int surfreactany;                 // 1 if any boundary has surf reactions

  int copy,copymode;            // 1 if copy of class (prevents deallocation of
                                //  base class when child copy is destroyed)

  int nregion;                      // # of defined Regions
  int maxregion;                    // max # regions can hold
  class Region **regions;           // list of defined Regions

  int surf_collide[6];              // index of SurfCollide model per face
  int surf_react[6];                // index of SurfReact model per face
                                    // for each bflag = SURFACE boundary

  Domain(class SPARTA *);
  virtual ~Domain();
  void init();
  void set_initial_box();
  void set_global_box();
  void set_boundary(int, char **);
  int periodic(int *);
  void boundary_modify(int, char **);
  virtual int collide(Particle::OnePart *&, int, int, double *, double &,
                      Particle::OnePart *&, int &);
  virtual void uncollide(int, double *);
  void add_region(int, char **);
  void delete_region(int, char **);
  int find_region(char *);
  void print_box(const char *);
};

}

#endif

/* ERROR/WARNING messages:

E: Axi-symmetry only allowed for 2d simulation

Self-explanatory.

E: Z dimension must be periodic for 2d simulation

Self-explanatory.

E: Box boundary not assigned a surf_collide ID

Any box boundary of style "s" must be assigned to a surface collision
model via the bound_modify command, before a simulation is performed.

E: Grid cutoff is longer than box length in a periodic dimension

This is not allowed.  Reduce the size of the cutoff specified by the
global gridcut command.

E: Box bounds are invalid

The box boundaries specified in the read_data file are invalid.  The
lo value must be less than the hi value for all 3 dimensions.

E: Boundary command after simulation box is defined

The boundary command cannot be used after a read_data, read_restart,
or create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Only ylo boundary can be axi-symmetric

Self-explanatory.  See the boundary doc page for more details.

E: Y cannot be periodic for axi-symmetric

Self-explanatory.  See the boundary doc page for more details.

E: Both sides of boundary must be periodic

Cannot specify a boundary as periodic only on the lo or hi side.  Must
be periodic on both sides.

E: Axi-symmetry is not yet supported in SPARTA

This error condition will be removed after axi-symmetry is
fully implemented.

E: Bound_modify surf requires wall be a surface

The box boundary must be of style "s" to be assigned a surface
collision model.

E: Bound_modify surf_collide ID is unknown

Self-explanatory.

E: Reuse of region ID

A region ID cannot be used twice.

E: Invalid region style

The choice of region style is unknown.

E: Delete region ID does not exist

Self-explanatory.

*/
