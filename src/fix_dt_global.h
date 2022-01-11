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

FixStyle(dt_global,FixDtGlobal)

#else

#ifndef SPARTA_FIX_DT_GLOBAL_H
#define SPARTA_FIX_DT_GLOBAL_H

#include "stdio.h"
#include "fix.h"
#include "compute.h"

namespace SPARTA_NS {

class FixDtGlobal : public Fix {
 public:
  FixDtGlobal(class SPARTA *, int, char **);
  ~FixDtGlobal();
  int setmask();
  void init();
  void end_of_step();
  virtual void reallocate();

protected:
  int nglocal;
  int lambdawhich,usqwhich,vsqwhihc,wsqwhich,tempwhich;
  char *id_lambda,*id_usq,*id_vsq,*id_wsq,*id_temp;
  int lambdaindex, usqindex,vsqindex,wsqindex,tempindex;
  class Compute *clambda;
  class Fix *flambda,*fusq,*fvsq,*fwsq,*ftemp;
  double *lambda,*usq,*vsq,*wsq,*temp;

 private:
  int me;
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
