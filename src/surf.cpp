/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#include "surf.h"
#include "memory.h"

using namespace DSMC_NS;

/* ---------------------------------------------------------------------- */

Surf::Surf(DSMC *dsmc) : Pointers(dsmc)
{
  surf_exist = 0;

  npoint = nline = ntri = 0;
  pts = NULL;
  lines = NULL;
  tris = NULL;
}

/* ---------------------------------------------------------------------- */

Surf::~Surf()
{
   memory->sfree(pts);
   memory->sfree(lines);
   memory->sfree(tris);
}

/* ---------------------------------------------------------------------- */

int Surf::add_id(char *idname)
{
  return 0;
}

/* ---------------------------------------------------------------------- */

void Surf::compute_line_normal(int nstart, int n)
{
}

/* ---------------------------------------------------------------------- */

void Surf::compute_tri_normal(int nstart, int n)
{
}

/* ---------------------------------------------------------------------- */

bigint Surf::memory_usage()
{
  bigint bytes = 0;
  bytes += (bigint) npoint * sizeof(Point);
  bytes += (bigint) nline * sizeof(Line);
  bytes += (bigint) ntri * sizeof(Tri);
  return bytes;
}
