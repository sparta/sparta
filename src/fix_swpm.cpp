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
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{INT,DOUBLE}; 

/* ---------------------------------------------------------------------- */

FixSWPM::FixSWPM(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg != 2) error->all(FLERR,"Illegal fix swpm command");

  flag_update_custom = 1;

  index_swpm = particle->find_custom((char *) "sweight");
  if (index_swpm > 0)
    error->all(FLERR,"Fix swpm custom attribute already exists");

  index_swpm = particle->add_custom((char *) "sweight",DOUBLE,0);

  sweight_new = atof(arg[1]);
  printf("%4.3e\n", sweight_new);
  error->one(FLERR,"ck");
}

/* ---------------------------------------------------------------------- */

FixSWPM::~FixSWPM()
{
  if (copy || copymode) return;

  particle->remove_custom(index_swpm);
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
  double *swpmweight = particle->edvec[particle->ewhich[index_swpm]];

  // conditional avoids weight reset during surface collisions
  // only new particles (zero weight) have their weight set as fnum

  if(swpmweight[index] == 0.0) swpmweight[index] = sweight_new;
}




