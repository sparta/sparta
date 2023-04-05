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

FixStyle(ave/histo/weight,FixAveHistoWeight)

#else

#ifndef SPARTA_FIX_AVE_HISTO_WEIGHT_H
#define SPARTA_FIX_AVE_HISTO_WEIGHT_H

#include <stdio.h>
#include "fix_ave_histo.h"

namespace SPARTA_NS {

class FixAveHistoWeight : public FixAveHisto {
 public:
  FixAveHistoWeight(class SPARTA *, int, char **);
  ~FixAveHistoWeight();

 private:
  int stridewt;

  // local method

  void bin_one_weight(double, double);

  // override these methods from ave/histo to use weights

  void bin_one(double);
  void bin_vector(int, double *, int);
  void bin_particles(int, int);
  void bin_particles(double *, int);
  void bin_grid_cells(double *, int);

  void calculate_weights();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Fix ave/histo/weight value and weight vector lengths do not match

Self-explanatory.

E: Invalid timestep reset for fix ave/histo

Resetting the timestep has invalidated the sequence of timesteps this
fix needs to process.

E: Error writing out histogram data

Something in the output to the file triggered an error.

*/
