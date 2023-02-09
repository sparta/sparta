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
#include "region_union.h"
#include "domain.h"
#include "error.h"

using namespace SPARTA_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

RegUnion::RegUnion(SPARTA *sparta, int narg, char **arg) :
  Region(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal region command");
  int n = atoi(arg[2]);
  if (n < 2) error->all(FLERR,"Illegal region command");
  options(narg-(n+3),&arg[n+3]);

  // build list of regions to union

  list = new int[n];
  nregion = 0;

  int iregion;
  for (int iarg = 0; iarg < n; iarg++) {
    iregion = domain->find_region(arg[iarg+3]);
    if (iregion == -1)
      error->all(FLERR,"Region union region ID does not exist");
    list[nregion++] = iregion;
  }

  // extent of union of regions
  // has bounding box if interior and all sub-regions have bounding box

  Region **regions = domain->regions;

  bboxflag = 1;
  for (int ilist = 0; ilist < nregion; ilist++)
    if (regions[list[ilist]]->bboxflag == 0) bboxflag = 0;
  if (!interior) bboxflag = 0;

  if (bboxflag) {
    extent_xlo = extent_ylo = extent_zlo = BIG;
    extent_xhi = extent_yhi = extent_zhi = -BIG;

    for (int ilist = 0; ilist < nregion; ilist++) {
      extent_xlo = MIN(extent_xlo,regions[list[ilist]]->extent_xlo);
      extent_ylo = MIN(extent_ylo,regions[list[ilist]]->extent_ylo);
      extent_zlo = MIN(extent_zlo,regions[list[ilist]]->extent_zlo);
      extent_xhi = MAX(extent_xhi,regions[list[ilist]]->extent_xhi);
      extent_yhi = MAX(extent_yhi,regions[list[ilist]]->extent_yhi);
      extent_zhi = MAX(extent_zhi,regions[list[ilist]]->extent_zhi);
    }
  }
}

/* ---------------------------------------------------------------------- */

RegUnion::~RegUnion()
{
  delete [] list;
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is match() with any sub-region
   else inside = 0
------------------------------------------------------------------------- */

int RegUnion::inside(double *x)
{
  int ilist;
  Region **regions = domain->regions;
  for (ilist = 0; ilist < nregion; ilist++)
    if (regions[list[ilist]]->match(x)) break;

  if (ilist == nregion) return 0;
  return 1;
}
