/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef REACT_CLASS

ReactStyle(qk,ReactQK)

#else

#ifndef SPARTA_REACT_QK_H
#define SPARTA_REACT_QK_H

#include "react.h"
#include "particle.h"

namespace SPARTA_NS {

class ReactQK : public React {
 public:
  ReactQK(class SPARTA *, int, char **);
  ~ReactQK() {}
  void init();
  int attempt(Particle::OnePart *, Particle::OnePart *, 
              double, double, double, double &, int &);

 private:
  struct Reaction {
    char *formula;
    int outcome;
    int ctrlmode;
    int nproducts;
    int number; // ??
    int product_C;
    int product_D;
    int product_E;
    int product_F;
    int product_G;
    int product_H;
    double e_forward;
    double e_backward;
    double vibn_disposal_param;
    double lambda_0;
    double lambda_t;
    double lambda_cutoff;
    double steric_factor;
  };

  Reaction **reactions;
  int nreactions;
};

}

#endif
#endif
