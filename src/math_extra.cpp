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

#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "math_extra.h"

namespace MathExtra {

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   1 = lower bound, Nmax = upper bound
   5 possibilities:
     (1) i = i to i, (2) * = 1 to Nmax,
     (3) i* = i to Nmax, (4) *j = 1 to j, (5) i*j = i to j
   return nlo,nhi
   return 0 if successful
   return 1 if numeric values are out of lower/upper bounds
------------------------------------------------------------------------- */

int bounds(char *str, int nmax, int &nlo, int &nhi)
{
  char *ptr = strchr(str,'*');

  if (ptr == NULL) {
    nlo = nhi = atoi(str);
  } else if (strlen(str) == 1) {
    nlo = 1;
    nhi = nmax;
  } else if (ptr == str) {
    nlo = 1;
    nhi = atoi(ptr+1);
  } else if (strlen(ptr+1) == 0) {
    nlo = atoi(str);
    nhi = nmax;
  } else {
    nlo = atoi(str);
    nhi = atoi(ptr+1);
  }

  if (nlo < 1 || nhi > nmax) return 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

}
