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

#ifdef COMPUTE_CLASS

ComputeStyle(boundary,ComputeBoundary)

#else

#ifndef SPARTA_BOUNDARY_SURF_H
#define SPARTA_BOUNDARY_SURF_H

#include "compute.h"
#include "particle.h"

namespace SPARTA_NS {

class ComputeBoundary : public Compute {
 public:
  ComputeBoundary(class SPARTA *, int, char **);
  ~ComputeBoundary();
  void init();
  void compute_array();
  void clear();
  void boundary_tally(int, int, double *, Particle::OnePart *);

 private:
  int imix,nvalue,ngroup,ntotal,nrow;
  int *which;

  double normflux[6];            // per face normalization factors
  double **myarray;              // local accumulator array
};

}

#endif
#endif
