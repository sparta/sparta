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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "react_tce_kokkos.h"
#include "particle.h"
#include "collide.h"
#include "random_knuth.h"
#include "error.h"

// DEBUG
#include "update.h"

using namespace SPARTA_NS;

enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};   // other files

/* ---------------------------------------------------------------------- */

ReactTCEKokkos::ReactTCEKokkos(SPARTA *sparta, int narg, char **arg) :
  ReactBirdKokkos(sparta, narg, arg) {}

/* ---------------------------------------------------------------------- */

void ReactTCEKokkos::init()
{
  if (!collide || (strcmp(collide->style,"vss") != 0 && strcmp(collide->style,"vss/kk") != 0))
    error->all(FLERR,"React tce can only be used with collide vss");

  ReactBirdKokkos::init();

  vibstyle = collide->vibstyle;
}

/* ---------------------------------------------------------------------- */

double ReactTCEKokkos::bird_Evib(int nmode, double Tvib,
                            double vibtemp[],
                            double Evib)
{

  // COMPUTES f FOR NEWTON'S SEARCH METHOD OUTLINED IN "newtonTvib".

  double f = -Evib;
  double kb = 1.38064852e-23;

  for (int i = 0; i < nmode; i++) {
    f += (((kb*vibtemp[i])/(exp(vibtemp[i]/Tvib)-1)));
  }

  return f;

}


/* ---------------------------------------------------------------------- */

double ReactTCEKokkos::bird_dEvib(int nmode, double Tvib, double vibtemp[])
{

  // COMPUTES df FOR NEWTON'S SEARCH METHOD

  double df = 0.0;
  double kb = 1.38064852e-23;

  for (int i = 0; i < nmode; i++) {
    df += ((pow(vibtemp[i],2)*kb*exp(vibtemp[i]/Tvib))/(pow(Tvib,2)*pow(exp(vibtemp[i]/Tvib)-1,2)));
  }

  return df;

}

/* ---------------------------------------------------------------------- */

double ReactTCEKokkos::newtonTvib(const int nmode, double Evib, const double vibTemp[],
               double Tvib0,
               double tol,
               int nmax)
{


  // FUNCTION FOR CONVERTING VIBRATIONAL ENERGY TO VIBRATIONAL TEMPERATURE
  // Computes Tvib assuming the vibrational energy levels occupy a simple harmonic oscillator (SHO)
  // spacing.
  // Search for Tvib begins at some initial value "Tvib0" until the search reaches a tolerance level "tol".


  double f;
  double df;
  double Tvib, Tvib_prev;
  double err;
  int i;

  // Uses Newton's method to solve for a vibrational temperature given a
  // distribution of vibrational energy levels.

  // f and df are computed for Newton's search
  f = bird_Evib(nmode,Tvib0,vibTemp,Evib);
  df = bird_dEvib(nmode,Tvib0,vibTemp);

  // Update guess for Tvib and compute error
  Tvib = Tvib0 - (f/df);
  err = fabs(Tvib-Tvib0);

  i=2;

  // Continue to search for Tvib until the error is greater than the tolerance:
  while((err >= tol) && (i <= nmax))
  {

    Tvib_prev = Tvib;

    f = bird_Evib(nmode,Tvib,vibTemp,Evib);
    df = bird_dEvib(nmode,Tvib,vibTemp);

    Tvib = Tvib_prev-(f/df);
    err = fabs(Tvib-Tvib_prev);

    i=i+1;

  }

  return Tvib;

}
