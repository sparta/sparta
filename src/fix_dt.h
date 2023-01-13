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

#include "fix.h"
#include "compute.h"

namespace SPARTA_NS {

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
  int mode;
  char *id_lambda,*id_usq,*id_vsq,*id_wsq,*id_temp;
  class Compute *clambda;
  class Fix *flambda,*fusq,*fvsq,*fwsq,*ftemp;
  double *lambda,*usq,*vsq,*wsq,*temp;
  double min_species_mass;
  double dt_global_weight;
  double transit_fraction;
  double collision_fraction;
  double dt_global_calculated;
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
