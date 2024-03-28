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
#include "stdlib.h"
#include "string.h"
#include "fix_swpm.h"
#include "update.h"
#include "particle.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{INT,DOUBLE};                      // several files

/* ---------------------------------------------------------------------- */

FixSWPM::FixSWPM(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix ambipolar command");

  flag_update_custom = 1;
  flag_surf_react = 1;

  swpmindex = particle->find_custom((char *) "sweight");

  if (swpmindex > 0)
    error->all(FLERR,"Fix ambipolar custom attribute already exists");

  swpmindex = particle->add_custom((char *) "sweight",DOUBLE,0);
}

/* ---------------------------------------------------------------------- */

FixSWPM::~FixSWPM()
{
  if (copy || copymode) return;

  particle->remove_custom(swpmindex);
}

/* ---------------------------------------------------------------------- */

int FixSWPM::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSWPM::init()
{
  return;
}

/* ----------------------------------------------------------------------
   called when a particle with index is created. ignores everything
   but index
------------------------------------------------------------------------- */

void FixSWPM::update_custom(int index, double,
                                double, double,
                                double*)
{
  double *swpmweight = particle->edvec[particle->ewhich[swpmindex]];
  swpmweight[index] = update->fnum;
}




