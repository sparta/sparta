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

#ifdef COMPUTE_CLASS

ComputeStyle(dt/grid,ComputeDtGrid)

#else

#ifndef SPARTA_COMPUTE_DT_GRID_H
#define SPARTA_COMPUTE_DT_GRID_H

#include "compute.h"

namespace SPARTA_NS {

class ComputeDtGrid : public Compute {
 public:
  ComputeDtGrid(class SPARTA *, int, char **);
  ~ComputeDtGrid();
  void init();
  virtual void compute_per_grid();
  virtual void reallocate();
  bigint memory_usage();

 protected:

  double transit_fraction,collision_fraction;
  double min_species_mass;

  int groupbit;
  int nglocal;
  int lambda_which,temp_which,usq_which,vsq_which,wsq_which;
  int lambda_index,temp_index,usq_index,vsq_index,wsq_index;;

  char *id_lambda,*id_temp,*id_usq,*id_vsq,*id_wsq;

  class Compute *clambda,*ctemp,*cusq,*cvsq,*cwsq;
  class Fix *flambda,*ftemp,*fusq,*fvsq,*fwsq;

  double *lambda,*temp,*usq,*vsq,*wsq;
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
