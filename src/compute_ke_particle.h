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

#ifdef COMPUTE_CLASS

ComputeStyle(ke/particle,ComputeKEParticle)

#else

#ifndef DSMC_COMPUTE_KE_PARTICLE_H
#define DSMC_COMPUTE_KE_PARTICLE_H

#include "compute.h"

namespace DSMC_NS {

class ComputeKEParticle : public Compute {
 public:
  ComputeKEParticle(class DSMC *, int, char **);
  ~ComputeKEParticle();
  void init();
  void compute_per_particle();
  double memory_usage();

 private:
  int nmax;
  double *ke;
};

}

#endif
#endif
