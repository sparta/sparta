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

#ifndef SPARTA_FIX_EMIT_H
#define SPARTA_FIX_EMIT_H

#include "fix.h"

namespace SPARTA_NS {

class FixEmit : public Fix {
 public:
  FixEmit(class SPARTA *, int, char **);
  virtual ~FixEmit();
  int setmask();
  virtual void init();
  void start_of_step();
  double compute_vector(int);

  virtual void grid_changed();

 protected:
  int perspecies;
  class Region *region;
  class RanKnuth *random;
  int nsingle,ntotal;

  int ntask;           // # of insert tasks in underlying child class

  int active_current;  // set to 0 if grid cell data struct changes
                       // triggers rebuild of active cell list in child classes

  virtual void create_task(int) = 0;
  virtual void perform_task() = 0;

  void create_tasks();
  double mol_inflow(double, double, double);
  int subsonic_temperature_check(int, double);
  void options(int, char **);
  virtual int option(int, char **);
};

}

#endif

/* ERROR/WARNING messages:

*/
