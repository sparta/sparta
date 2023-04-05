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

#include "stdlib.h"
#include "string.h"
#include "region_block.h"
#include "domain.h"
#include "error.h"

using namespace SPARTA_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

RegBlock::RegBlock(SPARTA *sparta, int narg, char **arg) :
  Region(sparta, narg, arg)
{
  options(narg-8,&arg[8]);

  if (strcmp(arg[2],"INF") == 0) xlo = -BIG;
  else xlo = atof(arg[2]);
  if (strcmp(arg[3],"INF") == 0) xhi = BIG;
  else xhi = atof(arg[3]);

  if (strcmp(arg[4],"INF") == 0) ylo = -BIG;
  else ylo = atof(arg[4]);
  if (strcmp(arg[5],"INF") == 0) yhi = BIG;
  else yhi = atof(arg[5]);

  if (strcmp(arg[6],"INF") == 0) zlo = -BIG;
  else zlo = atof(arg[6]);
  if (strcmp(arg[7],"INF") == 0) zhi = BIG;
  else zhi = atof(arg[7]);

  // error check

  if (xlo > xhi || ylo > yhi || zlo > zhi)
    error->all(FLERR,"Illegal region block command");

  // extent of block

  if (interior) {
    bboxflag = 1;
    extent_xlo = xlo;
    extent_xhi = xhi;
    extent_ylo = ylo;
    extent_yhi = yhi;
    extent_zlo = zlo;
    extent_zhi = zhi;
  } else bboxflag = 0;
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegBlock::inside(double *x)
{
  if (x[0] >= xlo && x[0] <= xhi && x[1] >= ylo && x[1] <= yhi &&
      x[2] >= zlo && x[2] <= zhi)
    return 1;
  return 0;
}
