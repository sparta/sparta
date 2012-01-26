/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(inflow,FixInflow)

#else

#ifndef DSMC_FIX_INFLOW_H
#define DSMC_FIX_INFLOW_H

#include "stdio.h"
#include "fix.h"

namespace DSMC_NS {

class FixInflow : public Fix {
 public:
  FixInflow(class DSMC *, int, char **);
  ~FixInflow();
  int setmask();
  void init();
  void start_of_step();
  double compute_vector(int);

 private:
  int nevery,imix;
  int faces[6];
  int np,nonce,ntotal;
  int npercell,nthresh;

  struct CellFace {
    double lo[3];
    double hi[3];
    double normal[3];
    double ntarget;
    double *ntargetsp;
  };

  CellFace *cellface;
  int ncf;

  class RanPark *random;

  double mol_inflow(int, double);
};

}

#endif
#endif
