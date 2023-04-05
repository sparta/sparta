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

#ifndef SPARTA_MODIFY_KOKKOS_H
#define SPARTA_MODIFY_KOKKOS_H

#include "modify.h"

namespace SPARTA_NS {

class ModifyKokkos : public Modify {
 public:
  ModifyKokkos(class SPARTA *);
  ~ModifyKokkos();
  void start_of_step();
  void end_of_step();

  int pack_grid_one(int, char *, int);
  int unpack_grid_one(int, char *);
  void copy_grid_one(int, int);
  void add_grid_one();
  void reset_grid_count(int);
  void grid_changed();

  void update_custom(int, double, double, double, double *);
  void gas_react(int);
  void surf_react(Particle::OnePart *, int &, int &);

 private:
  class ParticleKokkos* particle_kk;
  class GridKokkos* grid_kk;
};

}

#endif

/* ERROR/WARNING messages:

E: Fix command before simulation box is defined

The fix command cannot be used before a read_data, read_restart, or
create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Replacing a fix, but new style != old style

A fix ID can be used a 2nd time, but only if the style matches the
previous fix.  In this case it is assumed you with to reset a fix's
parameters.  This error may mean you are mistakenly re-using a fix ID
when you do not intend to.

E: Invalid fix style

The choice of fix style is unknown.

E: Could not find fix ID to delete

Self-explanatory.

E: Reuse of compute ID

A compute ID cannot be used twice.

E: Invalid compute style

Self-explanatory.

E: Could not find compute ID to delete

Self-explanatory.

*/
