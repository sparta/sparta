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
#include "region_intersect.h"
#include "domain.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

RegIntersect::RegIntersect(SPARTA *sparta, int narg, char **arg) :
  Region(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal region command");
  int n = atoi(arg[2]);
  if (n < 2) error->all(FLERR,"Illegal region command");
  options(narg-(n+3),&arg[n+3]);

  // build list of regions to intersect

  list = new int[n];
  nregion = 0;

  int iregion;
  for (int iarg = 0; iarg < n; iarg++) {
    iregion = domain->find_region(arg[iarg+3]);
    if (iregion == -1)
      error->all(FLERR,"Region intersect region ID does not exist");
    list[nregion++] = iregion;
  }

  // extent of intersection of regions
  // has bounding box if interior and any sub-region has bounding box

  Region **regions = domain->regions;

  bboxflag = 0;
  for (int ilist = 0; ilist < nregion; ilist++)
    if (regions[list[ilist]]->bboxflag == 1) bboxflag = 1;
  if (!interior) bboxflag = 0;

  if (bboxflag) {
    int first = 1;
    for (int ilist = 0; ilist < nregion; ilist++) {
      if (regions[list[ilist]]->bboxflag == 0) continue;
      if (first) {
        extent_xlo = regions[list[ilist]]->extent_xlo;
        extent_ylo = regions[list[ilist]]->extent_ylo;
        extent_zlo = regions[list[ilist]]->extent_zlo;
        extent_xhi = regions[list[ilist]]->extent_xhi;
        extent_yhi = regions[list[ilist]]->extent_yhi;
        extent_zhi = regions[list[ilist]]->extent_zhi;
        first = 0;
      }

      extent_xlo = MAX(extent_xlo,regions[list[ilist]]->extent_xlo);
      extent_ylo = MAX(extent_ylo,regions[list[ilist]]->extent_ylo);
      extent_zlo = MAX(extent_zlo,regions[list[ilist]]->extent_zlo);
      extent_xhi = MIN(extent_xhi,regions[list[ilist]]->extent_xhi);
      extent_yhi = MIN(extent_yhi,regions[list[ilist]]->extent_yhi);
      extent_zhi = MIN(extent_zhi,regions[list[ilist]]->extent_zhi);
    }
  }
}

/* ---------------------------------------------------------------------- */

RegIntersect::~RegIntersect()
{
  delete [] list;
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is match() with all sub-regions
   else inside = 0
------------------------------------------------------------------------- */

int RegIntersect::inside(double *x)
{
  int ilist;
  Region **regions = domain->regions;
  for (ilist = 0; ilist < nregion; ilist++)
    if (!regions[list[ilist]]->match(x)) break;

  if (ilist == nregion) return 1;
  return 0;
}
