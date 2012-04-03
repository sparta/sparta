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
#include "math_extra.h"
#include "memory.h"

using namespace DSMC_NS;
using namespace MathExtra;

/* ---------------------------------------------------------------------- */

Surf::Surf(DSMC *dsmc) : Pointers(dsmc)
{
  surf_exist = 0;

  nsurf = 0;
  ids = NULL;

  npoint = nline = ntri = 0;
  pts = NULL;
  lines = NULL;
  tris = NULL;
}

/* ---------------------------------------------------------------------- */

Surf::~Surf()
{
  for (int i = 0; i < nsurf; i++) delete [] ids[i];
  memory->sfree(ids);

  memory->sfree(pts);
  memory->sfree(lines);
  memory->sfree(tris);
}

/* ----------------------------------------------------------------------
   if idnew is already in ids list, return index
   else add it to ids list and return new index
------------------------------------------------------------------------- */

int Surf::add_id(char *idnew)
{
  for (int i = 0; i < nsurf; i++)
    if (strcmp(idnew,ids[i]) == 0) return i;

  ids = (char **) memory->srealloc(ids,(nsurf+1)*sizeof(char *),"surf:id");
  int n = strlen(idnew) + 1;
  ids[nsurf] = new char[n];
  strcpy(ids[nsurf],idnew);
  nsurf++;

  return nsurf-1;
}

/* ----------------------------------------------------------------------
   compute unit outward normal vectors of N lines starting at Nstart
   outward normal = +z axis x (p2-p1)
------------------------------------------------------------------------- */

void Surf::compute_line_normal(int nstart, int n)
{
  int p1,p2;
  double z[3],delta[3],norm[3];

  z[0] = 0.0; z[1] = 0.0; z[2] = 1.0;

  int m = nstart;
  for (int i = 0; i < n; i++) {
    p1 = lines[m].p1;
    p2 = lines[m].p2;
    MathExtra::sub3(pts[p2].x,pts[p1].x,delta);
    MathExtra::cross3(z,delta,norm);
    MathExtra::norm3(norm);
    lines[m].norm[0] = norm[0];
    lines[m].norm[1] = norm[1];
    m++;
  }
}

/* ----------------------------------------------------------------------
   compute unit outward normal vectors of N triangles starting at Nstart
   outward normal = (p2-p1) x (p3-p1)
------------------------------------------------------------------------- */

void Surf::compute_tri_normal(int nstart, int n)
{
  int p1,p2,p3;
  double delta12[3],delta13[3];

  int m = nstart;
  for (int i = 0; i < n; i++) {
    p1 = tris[m].p1;
    p2 = tris[m].p2;
    p3 = tris[m].p3;
    MathExtra::sub3(pts[p2].x,pts[p1].x,delta12);
    MathExtra::sub3(pts[p3].x,pts[p1].x,delta13);
    MathExtra::cross3(delta12,delta13,tris[m].norm);
    MathExtra::norm3(tris[m].norm);
    m++;
  }
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
