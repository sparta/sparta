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

#include "spatype.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "math_extra.h"

using namespace SPARTA_NS;

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

/* ----------------------------------------------------------------------
   convert a 64-bit integer into a string like "9.36B"
   K = thousand, M = million, B = billion, T = trillion, P = peta, E = exa
   for easier-to-understand output
------------------------------------------------------------------------- */

char *num2str(bigint n, char *outstr)
{
  if (n < 100000) sprintf(outstr,"(%1.3gK)",1.0e-3*n);
  else if (n < 1000000000) sprintf(outstr,"(%1.3gM)",1.0e-6*n);
  else if (n < 1000000000000) sprintf(outstr,"(%1.3gB)",1.0e-9*n);
  else if (n < 1000000000000000) sprintf(outstr,"(%1.3gT)",1.0e-12*n);
  else if (n < 1000000000000000000) sprintf(outstr,"(%1.3gP)",1.0e-15*n);
  else sprintf(outstr,"(%1.3gE)",1.0e-18*n);
  return outstr;
}

}
