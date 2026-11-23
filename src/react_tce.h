/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef REACT_CLASS

ReactStyle(tce,ReactTCE)

#else

#ifndef SPARTA_REACT_TCE_H
#define SPARTA_REACT_TCE_H

#include "react_bird.h"
#include "particle.h"

namespace SPARTA_NS {

class ReactTCE : public ReactBird {
 public:
  ReactTCE(class SPARTA *, int, char **);
  void init();
  void end_of_run();
  int attempt(Particle::OnePart *, Particle::OnePart *,
              double, double, double, double &, int &);

  double newtonTvib(int nmode, double Evib,
                      double VibTemp[],
                      double Tvib0,
                      double tol,
                      int nmax);

  double bird_dEvib(int nmode, double Tvib,
                  double VibTemp[]);

  double bird_Evib(int nmode, double Tvib,
                 double VibTemp[],
                 double Evib);

 protected:
  int prob_neg_flag;       // 1 if a react_prob < 0 occurred on this proc
                           //   during the current run
  int prob_big_index;      // rlist index of the first reaction whose summed
                           //   react_prob exceeded 1 on this proc during
                           //   the current run, -1 if none
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

W: ... TCE reaction(s) have a temperature exponent below the bound
-(z+3/2) for which the TCE model can reproduce their Arrhenius rate ...

For each listed reaction the Arrhenius rate cannot be represented by
the TCE model for the energy and vibrational models in use: the
argument of the gamma function in the TCE reaction probability is
non-positive, so the gamma function is clamped and the reaction rate
is incorrect.  With discrete vibration this only happens when the
vibrational energy is small.  Use reaction coefficients whose
temperature exponent is above the printed bound, e.g. by increasing
the number of internal degrees of freedom z where that is physically
valid.  This warning is printed once, on the first run.

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
