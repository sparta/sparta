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
#include "react.h"
#include "update.h"
#include "comm.h"
#include "input.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

React::React(SPARTA *sparta, int, char **arg) : Pointers(sparta)
{
  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  recombflag_user = 1;
  recomb_boost = 1000.0;
  recomb_boost_inverse = 0.001;

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);

  copy = copymode = 0;
}

/* ---------------------------------------------------------------------- */

React::~React()
{
  if (copy) return;

  delete [] style;
  delete random;
}

/* ---------------------------------------------------------------------- */

void React::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal react_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"recomb") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal react_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) recombflag_user = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) recombflag_user = 0;
      else error->all(FLERR,"Illegal react_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"rboost") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal react_modify command");
      recomb_boost = input->numeric(FLERR,arg[iarg+1]);
      if (recomb_boost < 1.0) error->all(FLERR,"Illegal react_modify command");
      recomb_boost_inverse = 1.0 / recomb_boost;
      iarg += 2;
    } else error->all(FLERR,"Illegal react_modify command");
  }
}
