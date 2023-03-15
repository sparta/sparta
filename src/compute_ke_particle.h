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

ComputeStyle(ke/particle,ComputeKEParticle)

#else

#ifndef SPARTA_COMPUTE_KE_PARTICLE_H
#define SPARTA_COMPUTE_KE_PARTICLE_H

#include "compute.h"

namespace SPARTA_NS {

class ComputeKEParticle : public Compute {
 public:
  ComputeKEParticle(class SPARTA *, int, char **);
  ~ComputeKEParticle();
  void init();
  void compute_per_particle();
  bigint memory_usage();

 protected:
  int nmax;

 private:
  double *ke;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

W: More than one compute ke/particle

This may be inefficient since each such compute stores a vector
of length equal to the number of particles.

*/
