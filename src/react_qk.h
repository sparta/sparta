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

#ifdef REACT_CLASS

ReactStyle(qk,ReactQK)

#else

#ifndef SPARTA_REACT_QK_H
#define SPARTA_REACT_QK_H

#include "react_bird.h"
#include "particle.h"

namespace SPARTA_NS {

class ReactQK : public ReactBird {
 public:
  ReactQK(class SPARTA *, int, char **);
  void init();
  int attempt(Particle::OnePart *, Particle::OnePart *,
              double, double, double, double &, int &);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: React tce can only be used with collide vss

Self-explanatory.

E: Ionization and recombination reactions are not yet implemented

This error conditions will be removed after those reaction styles are
fully implemented.

E: Unknown outcome in reaction

The specified type of the reaction is not encoded in the reaction
style.

E: Cannot open reaction file %s

Self-explanatory.

E: Invalid reaction formula in file

Self-explanatory.

E: Invalid reaction type in file

Self-explanatory.

E: Invalid reaction style in file

Self-explanatory.

E: Invalid reaction coefficients in file

Self-explanatory.

*/
