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

#ifdef COMPUTE_CLASS

ComputeStyle(surf,ComputeSurf)

#else

#ifndef SPARTA_COMPUTE_SURF_H
#define SPARTA_COMPUTE_SURF_H

#include "compute.h"

namespace SPARTA_NS {

class ComputeSurf : public Compute {
 public:
  ComputeSurf(class SPARTA *, int, char **);
  ~ComputeSurf();
  void init();
  void compute_per_surf();
  void clear();
  void tally(int isurf, int ispecies, double *v);
  double *normptr(int);
  bigint memory_usage();

 private:
  int imix,nvalue,ngroup,ntotal;
  int *which;



  int nbounce;






  int **value_norm_style;       // I,J = norm style of Jth value in Ith group
  double **norm_count;          // per-group ptr to norm vector, by count
  double **norm_mass;           // per-group ptr to norm vector, by mass
  double **norm_temp;           // per-group ptr to norm vector, by temperature
};

}

#endif
#endif
