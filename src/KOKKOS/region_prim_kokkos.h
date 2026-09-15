/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_REGION_PRIM_KOKKOS_H
#define SPARTA_REGION_PRIM_KOKKOS_H

#include "kokkos_type.h"

namespace SPARTA_NS {

// Region::inside() is a host virtual, and virtual dispatch is not available
//   inside a device kernel.  Rather than have every consumer carry a KKCopy
//   of each concrete region type and switch on a style string -- which is
//   what update/emit used to do, and what capped the number of regions a
//   run could use -- each Kokkos region flattens itself into a small
//   device-resident POSTFIX (RPN) token stream, which a kernel walks with a
//   fixed-depth boolean stack and no dispatch at all.
//
// a token is either
//   RKK_TOK_PRIM       push region_prim_match_kk(token.prim) onto the stack
//   RKK_TOK_UNION      pop 2, push (a || b)
//   RKK_TOK_INTERSECT  pop 2, push (a && b)
//   RKK_TOK_NOT        negate the top of stack
//
// a primitive region flattens to a single RKK_TOK_PRIM token whose own
//   interior/exterior sense is already folded into the primitive (the
//   !(inside ^ interior) of Region::match()).
// region union / region intersect flatten to the concatenation of their
//   sub-regions' streams, left-folded with one op token after each
//   sub-region past the first, followed by one RKK_TOK_NOT token when the
//   composite itself is an exterior region (interior == 0), since
//   !(hit ^ 0) == !hit and !(hit ^ 1) == hit.  because a sub-region
//   contributes a whole sub-stream rather than a single entry, composites
//   nest to arbitrary depth as long as the stack bound below holds.
//
// an op token carries an unused RegionPrimKK payload.  that wastes a little
//   memory per op, but keeps the whole program in ONE view, so a consumer
//   still holds a single DualView plus a token count.

enum{RKK_BLOCK,RKK_CYLINDER,RKK_PLANE,RKK_SPHERE};
enum{RKK_TOK_PRIM,RKK_TOK_UNION,RKK_TOK_INTERSECT,RKK_TOK_NOT};

// max boolean stack depth a token stream may require.  a stream that needs
//   more than this is rejected on the host at flatten time (see
//   region_token_depth() below and its callers) -- never truncated.
// a flat composite of any number of primitives needs depth 2; depth D
//   allows e.g. D-1 levels of right-nested composites.  16 is far past any
//   region tree an input script is likely to build, and costs the kernel
//   only a 16-int (64 byte) per-thread scratch array on the paths that are
//   not the single-primitive fast path.

enum{RKK_MAX_STACK = 16};

struct RegionPrimKK {
  int style;              // one of RKK_*
  int interior;           // this sub-region's own interior/exterior sense
  int axis;               // cylinder only: 0/1/2 for x/y/z

  // style-specific parameters, packed:
  //   BLOCK     a..f = xlo,xhi,ylo,yhi,zlo,zhi
  //   CYLINDER  a,b = c1,c2   c = radius   d,e = lo,hi
  //   PLANE     a,b,c = point   n0,n1,n2 = normal
  //   SPHERE    a,b,c = center  d = radius

  double a,b,c,d,e,f;
  double n0,n1,n2;
};

struct RegionTokenKK {
  int type;               // one of RKK_TOK_*
  RegionPrimKK prim;      // meaningful only when type == RKK_TOK_PRIM
};

typedef Kokkos::DualView<RegionTokenKK*, DeviceType::array_layout, DeviceType>
  tdual_region_token_1d;
typedef tdual_region_token_1d::t_dev t_region_token_1d;

/* ----------------------------------------------------------------------
   does x,y,z match a single flattened sub-region
   mirrors Region::match(): !(inside ^ interior)
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
int region_prim_match_kk(const RegionPrimKK &p,
                         const double x, const double y, const double z)
{
  int inside = 0;

  if (p.style == RKK_BLOCK) {
    if (x >= p.a && x <= p.b && y >= p.c && y <= p.d && z >= p.e && z <= p.f)
      inside = 1;

  } else if (p.style == RKK_CYLINDER) {
    double del1,del2;
    if (p.axis == 0) { del1 = y - p.a; del2 = z - p.b; }
    else if (p.axis == 1) { del1 = x - p.a; del2 = z - p.b; }
    else { del1 = x - p.a; del2 = y - p.b; }
    const double dist = sqrt(del1*del1 + del2*del2);
    const double along = (p.axis == 0) ? x : ((p.axis == 1) ? y : z);
    if (dist <= p.c && along >= p.d && along <= p.e) inside = 1;

  } else if (p.style == RKK_PLANE) {
    const double dot = (x-p.a)*p.n0 + (y-p.b)*p.n1 + (z-p.c)*p.n2;
    if (dot >= 0.0) inside = 1;

  } else {   // RKK_SPHERE
    const double delx = x - p.a;
    const double dely = y - p.b;
    const double delz = z - p.c;
    if (sqrt(delx*delx + dely*dely + delz*delz) <= p.d) inside = 1;
  }

  return !(inside ^ p.interior);
}

/* ----------------------------------------------------------------------
   does x,y,z match a flattened region, i.e. evaluate its postfix token
     stream against a small boolean stack
   ntoken == 1 is the single-primitive fast path: no stack at all, exactly
     the work the flat representation used to do
   the stream is built and validated on the host (region_token_depth), so
     the stack can never overflow or underflow here; the guards below are
     belt-and-braces against a caller passing a stale ntoken
------------------------------------------------------------------------- */

template<class ViewType>
KOKKOS_INLINE_FUNCTION
int region_match_kk(const ViewType &d_tokens, const int ntoken,
                    const double x, const double y, const double z)
{
  if (ntoken == 1) return region_prim_match_kk(d_tokens[0].prim,x,y,z);
  if (ntoken <= 0) return 0;

  int stack[RKK_MAX_STACK];
  int nstack = 0;

  for (int i = 0; i < ntoken; i++) {
    const int type = d_tokens[i].type;

    if (type == RKK_TOK_PRIM) {
      if (nstack == RKK_MAX_STACK) return 0;
      stack[nstack++] = region_prim_match_kk(d_tokens[i].prim,x,y,z);

    } else if (type == RKK_TOK_NOT) {
      if (nstack < 1) return 0;
      stack[nstack-1] = !stack[nstack-1];

    } else {
      if (nstack < 2) return 0;
      const int b = stack[--nstack];
      const int a = stack[nstack-1];
      if (type == RKK_TOK_UNION) stack[nstack-1] = (a || b);
      else stack[nstack-1] = (a && b);
    }
  }

  return stack[0];
}

/* ----------------------------------------------------------------------
   host-side helpers used while a composite region builds its token stream
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   make sure k_tokens can hold n tokens, preserving the tokens already in it
   grows geometrically; the host view is the authoritative copy while a
     stream is being built, so only host data is carried over
------------------------------------------------------------------------- */

inline void region_token_grow(tdual_region_token_1d &k_tokens, const int n)
{
  const int nmax = (int) k_tokens.extent(0);
  if (nmax >= n) return;

  int nnew = nmax ? nmax : 8;
  while (nnew < n) nnew *= 2;

  tdual_region_token_1d k_new("region:tokens",nnew);
  for (int i = 0; i < nmax; i++)
    k_new.view_host()[i] = k_tokens.view_host()[i];
  k_tokens = k_new;
}

/* ----------------------------------------------------------------------
   boolean stack depth the first ntoken tokens of k_tokens will require
   returns -1 if the stream is malformed (stack underflow, or anything
     other than exactly one value left at the end)
   host only: called at flatten time so a too-deep or ill-formed region tree
     is an error, not a wrong answer in a kernel
------------------------------------------------------------------------- */

inline int region_token_depth(tdual_region_token_1d &k_tokens, const int ntoken)
{
  if (ntoken <= 0) return -1;
  if ((int) k_tokens.extent(0) < ntoken) return -1;

  int nstack = 0;
  int maxdepth = 0;

  for (int i = 0; i < ntoken; i++) {
    const int type = k_tokens.view_host()[i].type;

    if (type == RKK_TOK_PRIM) {
      nstack++;
      if (nstack > maxdepth) maxdepth = nstack;
    } else if (type == RKK_TOK_NOT) {
      if (nstack < 1) return -1;
    } else {
      if (nstack < 2) return -1;
      nstack--;
    }
  }

  if (nstack != 1) return -1;
  return maxdepth;
}

}

#endif
