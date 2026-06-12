/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   File adapted from LAMMPS (https://www.lammps.org), October 2024
   Ported to SPARTA by: Stan Moore (SNL)
   Original Author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(halt,FixHalt)

#else

#ifndef SPARTA_FIX_HALT_H
#define SPARTA_FIX_HALT_H

#include "fix.h"

namespace SPARTA_NS {

class FixHalt : public Fix {
 public:
  FixHalt(class SPARTA *, int, char **);
  ~FixHalt() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;
  void post_run() override;

 private:
  int attribute, operation, eflag, msgflag, ivar;
  bigint nextstep, thisstep;
  double value, tratio;
  char *idvar;

  double tlimit();
};

}    // namespace SPARTA_NS

#endif
#endif
