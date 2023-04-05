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
#include "random_knuth.h"
#include "stdlib.h"

using namespace SPARTA_NS;

#define IM 2147483647
#define MBIG 1000000000
#define MSEED 161803398
#define FAC (1.0/MBIG)

/* ----------------------------------------------------------------------
   Knuth RNG
   assume iseed is a positive int
------------------------------------------------------------------------ */

RanKnuth::RanKnuth(int iseed)
{
  seed = iseed;
  save = 0;
  initflag = 0;
}

/* ----------------------------------------------------------------------
   set seed to positive int
   assume 0.0 <= rseed < 1.0
------------------------------------------------------------------------ */

RanKnuth::RanKnuth(double rseed)
{
  seed = static_cast<int> (rseed*IM);
  if (seed == 0) seed = 1;
  save = 0;
  initflag = 0;
}

/* ----------------------------------------------------------------------
   reset seed to a positive int based on rseed and offset
   assume 0.0 <= rseed < 1.0 and offset is an int >= 0
   fmod() insures no overflow when static cast to int
   warmup the new RNG if requested
   typically used to setup one RN generator per proc or site or particle
------------------------------------------------------------------------ */

void RanKnuth::reset(double rseed, int offset, int warmup)
{
  seed = static_cast<int> (fmod(rseed*IM+offset,IM));
  if (seed < 0) seed = -seed;
  if (seed == 0) seed = 1;
  initflag = 0;
  for (int i = 0; i < warmup; i++) uniform();
}

/* ----------------------------------------------------------------------
   uniform RN
------------------------------------------------------------------------- */

double RanKnuth::uniform()
{
  int i,ii,k,mj,mk;
  double rn;

  if (!initflag) {
    initflag = 1;
    mj = labs(MSEED-labs(seed));
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    for (i = 1; i <= 54; i++) {
      ii = (21*i) % 55;
      ma[ii] = mk;
      mk = mj-mk;
      if (mk < 0) mk += MBIG;
      mj = ma[ii];
    }
    for (k=0; k<4; k++)
      for (i=1; i<=55; i++) {
        ma[i] -= ma[1+(i+30) % 55];
        if (ma[i] < 0) ma[i] += MBIG;
    }
    inext = 0;
    inextp = 31;
  }

  while (1) {
    if (++inext == 56) inext = 1;
    if (++inextp == 56) inextp = 1;
    mj = ma[inext] - ma[inextp];
    if (mj < 0) mj += MBIG;
    ma[inext] = mj;
    rn = mj*FAC;

    // make sure the random number is valid

    if (rn > 0.0 && rn < 1.0) break;
  }

  return rn;
}

/* ----------------------------------------------------------------------
   gaussian RN with zero mean and unit variance
------------------------------------------------------------------------- */

double RanKnuth::gaussian()
{
  double first,v1,v2,rsq,fac;

  if (!save) {
    while (1) {
      v1 = 2.0*uniform()-1.0;
      v2 = 2.0*uniform()-1.0;
      rsq = v1*v1 + v2*v2;
      if (rsq < 1.0 && rsq != 0.0) break;
    }
    fac = sqrt(-2.0*log(rsq)/rsq);
    second = v1*fac;
    first = v2*fac;
    save = 1;
  } else {
    first = second;
    save = 0;
  }
  return first;
}
