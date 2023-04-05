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

#ifdef FIX_CLASS

FixStyle(balance,FixBalance)

#else

#ifndef SPARTA_FIX_BALANCE_H
#define SPARTA_FIX_BALANCE_H

#include "fix.h"

namespace SPARTA_NS {

class FixBalance : public Fix {
 public:
  FixBalance(class SPARTA *, int, char **);
  ~FixBalance();
  int setmask();
  void init();
  virtual void end_of_step();
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

 private:
  int me,nprocs;
  double thresh;
  int bstyle,rcbwt,rcbflip;
  char eligible[4];
  double last,my_timer_cost;

  double imbnow;                // current imbalance factor
  double imbprev;               // imbalance factor before last rebalancing
  double imbfinal;              // imbalance factor after last rebalancing
  double maxperproc;            // max atoms or CPU cost on any processor

  class RanKnuth *random;
  class RCB *rcb;

  double imbalance_factor(double &);
  void timer_cost();
  void timer_cell_weights(double *&);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Cannot use non-rcb fix balance with a grid cutoff

This is because the load-balancing will generate a partitioning
of cells to processors that is dispersed and which will not work
with a grid cutoff >= 0.0.

*/
