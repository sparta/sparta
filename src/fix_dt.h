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

#ifdef FIX_CLASS

FixStyle(dt,FixDt)

#else

#ifndef SPARTA_FIX_DT_H
#define SPARTA_FIX_DT_H

#include "stdio.h"
#include "fix.h"
#include "compute.h"

namespace SPARTA_NS {

enum class FIXMODE : int {NONE, WARN, USE_CALCULATED_GLOBAL_DT, USE_CALCULATED_CELL_DT};

class FixDt : public Fix {
 public:
  FixDt(class SPARTA *, int, char **);
  ~FixDt();
  int setmask();
  void init();
  void end_of_step();
  double compute_scalar();
  virtual void reallocate();

 protected:
  int me;
  int nglocal;
  int lambdawhich,usqwhich,vsqwhihc,wsqwhich,tempwhich;
  int lambdaindex, usqindex,vsqindex,wsqindex,tempindex;
  int tvar,txvar,tyvar,tzvar,tempflag;
  int vxvar,vyvar,vzvar,vvarx,vvary,vvarz,velflag;
  int imix;
  char *id_lambda,*id_usq,*id_vsq,*id_wsq,*id_temp;
  char *tstr,*txstr,*tystr,*tzstr;
  char *vxstr,*vystr,*vzstr,*vstrx,*vstry,*vstrz;
  class Compute *clambda;
  class Fix *flambda,*fusq,*fvsq,*fwsq,*ftemp;
  double *lambda,*usq,*vsq,*wsq,*temp;
  double min_species_mass;
  double dt_global_weight;
  double dt_global_calculated=0.;
  FIXMODE mode = FIXMODE::NONE;

  double temperature_variable(double *);
  void velocity_variable(double *, double *, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/
