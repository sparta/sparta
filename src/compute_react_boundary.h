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

ComputeStyle(react/boundary,ComputeReactBoundary)

#else

#ifndef SPARTA_REACT_BOUNDARY_SURF_H
#define SPARTA_REACT_BOUNDARY_SURF_H

#include "compute.h"
#include "particle.h"

namespace SPARTA_NS {

class ComputeReactBoundary : public Compute {
 public:
  ComputeReactBoundary(class SPARTA *, int, char **);
  ~ComputeReactBoundary();
  virtual void init();
  virtual void compute_array();
  virtual void clear();
  virtual void boundary_tally(int, int, int, Particle::OnePart *,
                              Particle::OnePart *, Particle::OnePart *);

 protected:
  int isr,ntotal,nrow,rpflag;
  int *surf_react;

  int **reaction2col;      // 1 if ireaction triggers tally for icol

  double **myarray;        // local accumulator array
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
