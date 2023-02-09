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
  ComputeBoundary(class SPARTA* sparta) : Compute(sparta) {}
  ~ComputeBoundary();
  virtual void init();
  virtual void compute_array();
  virtual void clear();
  virtual void boundary_tally(int, int, int, Particle::OnePart *,
                              Particle::OnePart *, Particle::OnePart *);

 protected:
  int imix,nvalue,ngroup,ntotal,nrow;
  int *which;

  int weightflag;                // 1 if cell weighting is enabled
  double weight;                 // particle weight, based on initial cell

  double normflux[6];            // per face normalization factors
  double **myarray;              // local accumulator array
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute boundary mixture ID does not exist

Self-explanatory.

E: Number of groups in compute boundary mixture has changed

This mixture property cannot be changed after this compute command is
issued.

*/
