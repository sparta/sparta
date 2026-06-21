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
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
FixStyle(controller,FixController)
#else

#ifndef SPARTA_FIX_CONTROLLER_H
#define SPARTA_FIX_CONTROLLER_H

#include "fix.h"

namespace SPARTA_NS {

class FixController : public Fix {
 public:
  FixController(class SPARTA *, int, char **);
  ~FixController() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;
  double compute_vector(int) override;

 private:
  double kp, ki, kd, alpha, tau;
  double setpoint;
  int pvwhich, pvindex;
  char *pvID, *cvID;
  int firsttime;

  double control, err, olderr, deltaerr, sumerr;

  class Compute *pcompute;
  class Fix *pfix;
  int pvar, cvar;
};

}    // namespace SPARTA_NS

#endif
#endif
