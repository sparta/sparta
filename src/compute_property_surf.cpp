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

#include "string.h"
#include "compute_property_surf.h"
#include "surf.h"
#include "update.h"
#include "domain.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputePropertySurf::ComputePropertySurf(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute property/surf command");

  dimension = domain->dimension;
  
  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute surf group ID does not exist");
  groupbit = surf->bitmask[igroup];

  nvalues = narg - 3;

  // parse input values
  // customize a new keyword by adding to if statement

  pack_choice = new FnPtrPack[nvalues];

  int i;
  for (int iarg = 3; iarg < narg; iarg++) {
    i = iarg-3;

    // check for invalid fields in 2d
    
    if (dimension == 2)
      if ((strcmp(arg[iarg],"v1z") == 0) || (strcmp(arg[iarg],"v2z") == 0) ||
	  (strcmp(arg[iarg],"v3x") == 0) || (strcmp(arg[iarg],"v3y") == 0) ||
	  (strcmp(arg[iarg],"v3z") == 0) || (strcmp(arg[iarg],"zc") == 0) ||
	  (strcmp(arg[iarg],"normz") == 0))
        error->all(FLERR,"Invalid compute property/surf field for 2d simulation");
    
    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_id;

    } else if (strcmp(arg[iarg],"v1x") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_v1x;
    } else if (strcmp(arg[iarg],"v1y") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_v1y;
    } else if (strcmp(arg[iarg],"v1z") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_v1z;
    } else if (strcmp(arg[iarg],"v2x") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_v2x;
    } else if (strcmp(arg[iarg],"v2y") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_v2y;
    } else if (strcmp(arg[iarg],"v2z") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_v2z;
    } else if (strcmp(arg[iarg],"v3x") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_v3x;
    } else if (strcmp(arg[iarg],"v3y") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_v3y;
    } else if (strcmp(arg[iarg],"v3z") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_v3z;

    } else if (strcmp(arg[iarg],"xc") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_xc;
    } else if (strcmp(arg[iarg],"yc") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_yc;
    } else if (strcmp(arg[iarg],"zc") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_zc;
      
    } else if (strcmp(arg[iarg],"area") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_area;
    } else if (strcmp(arg[iarg],"normx") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_normx;
    } else if (strcmp(arg[iarg],"normy") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_normy;
    } else if (strcmp(arg[iarg],"normz") == 0) {
      pack_choice[i] = &ComputePropertySurf::pack_normz;

    } else error->all(FLERR,"Invalid keyword in compute property/surf command");
  }

  
  per_surf_flag = 1;
  if (nvalues == 1) size_per_surf_cols = 0;
  else size_per_surf_cols = nvalues;

  nslocal = 0;
  vector_surf = NULL;
  array_surf = NULL;
}

/* ---------------------------------------------------------------------- */

ComputePropertySurf::~ComputePropertySurf()
{
  if (copymode) return;
  delete [] pack_choice;
  memory->destroy(vector_surf);
  memory->destroy(array_surf);
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::init()
{
  if (!surf->exist)
    error->all(FLERR,"Cannot use compute property/surf when surfs do not exist");
  if (surf->implicit)
    error->all(FLERR,"Cannot use compute property/surf with implicit surfs");

  distributed = surf->distributed;

  if (surf->nlocal == nslocal) return;

  nslocal = surf->nlocal;
  if (nvalues == 1) {
    memory->destroy(vector_surf);
    memory->create(vector_surf,nslocal,"property/surf:vector_surf");
  } else {
    memory->destroy(array_surf);
    memory->create(array_surf,nslocal,nvalues,"property/surf:array_surf");
  }

}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::compute_per_surf()
{
  invoked_per_surf = update->ntimestep;

  // fill vector or array with per-surf values

  if (nvalues == 1) {
    buf = vector_surf;
    (this->*pack_choice[0])(0);
  } else {
    if (nslocal) buf = &array_surf[0][0];
    else buf = NULL;
    for (int n = 0; n < nvalues; n++)
      (this->*pack_choice[n])(n);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local surf-based array
------------------------------------------------------------------------- */

bigint ComputePropertySurf::memory_usage()
{
  bigint bytes;
  bytes = nvalues*nslocal * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword compute property/surf can output
   the surf property is packed into buf starting at n with stride nvalues
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

void ComputePropertySurf::pack_id(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nslocal; i++) {
      if (lines[i].mask & groupbit) buf[n] = lines[i].id;
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nslocal; i++) {
      if (tris[i].mask & groupbit) buf[n] = tris[i].id;
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_v1x(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nslocal; i++) {
      if (lines[i].mask & groupbit) buf[n] = lines[i].p1[0];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nslocal; i++) {
      if (tris[i].mask & groupbit) buf[n] = tris[i].p1[0];
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_v1y(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nslocal; i++) {
      if (lines[i].mask & groupbit) buf[n] = lines[i].p1[1];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nslocal; i++) {
      if (tris[i].mask & groupbit) buf[n] = tris[i].p1[1];
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_v1z(int n)
{
  Surf::Tri *tris;
  if (distributed) tris = surf->mytris;
  else tris = surf->tris;
  for (int i = 0; i < nslocal; i++) {
    if (tris[i].mask & groupbit) buf[n] = tris[i].p1[2];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_v2x(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nslocal; i++) {
      if (lines[i].mask & groupbit) buf[n] = lines[i].p2[0];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nslocal; i++) {
      if (tris[i].mask & groupbit) buf[n] = tris[i].p2[0];
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_v2y(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nslocal; i++) {
      if (lines[i].mask & groupbit) buf[n] = lines[i].p2[1];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nslocal; i++) {
      if (tris[i].mask & groupbit) buf[n] = tris[i].p2[1];
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_v2z(int n)
{
  Surf::Tri *tris;
  if (distributed) tris = surf->mytris;
  else tris = surf->tris;
  for (int i = 0; i < nslocal; i++) {
    if (tris[i].mask & groupbit) buf[n] = tris[i].p2[2];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_v3x(int n)
{
  Surf::Tri *tris;
  if (distributed) tris = surf->mytris;
  else tris = surf->tris;
  for (int i = 0; i < nslocal; i++) {
    if (tris[i].mask & groupbit) buf[n] = tris[i].p3[0];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_v3y(int n)
{
  Surf::Tri *tris;
  if (distributed) tris = surf->mytris;
  else tris = surf->tris;
  for (int i = 0; i < nslocal; i++) {
    if (tris[i].mask & groupbit) buf[n] = tris[i].p1[1];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_v3z(int n)
{
  Surf::Tri *tris;
  if (distributed) tris = surf->mytris;
  else tris = surf->tris;
  for (int i = 0; i < nslocal; i++) {
    if (tris[i].mask & groupbit) buf[n] = tris[i].p1[2];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_xc(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nslocal; i++) {
      if (lines[i].mask & groupbit)
	buf[n] = 0.5 * (lines[i].p1[0] + lines[i].p2[0]);
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nslocal; i++) {
      if (tris[i].mask & groupbit)
	buf[n] = THIRD * (tris[i].p1[0] + tris[i].p2[0] + tris[i].p3[0]);
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_yc(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nslocal; i++) {
      if (lines[i].mask & groupbit)
	buf[n] = 0.5 * (lines[i].p1[1] + lines[i].p2[1]);
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nslocal; i++) {
      if (tris[i].mask & groupbit)
	buf[n] = THIRD * (tris[i].p1[1] + tris[i].p2[1] + tris[i].p3[1]);
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_zc(int n)
{
  Surf::Tri *tris;
  if (distributed) tris = surf->mytris;
  else tris = surf->tris;
  for (int i = 0; i < nslocal; i++) {
    if (tris[i].mask & groupbit)
      buf[n] = THIRD * (tris[i].p1[2] + tris[i].p2[2] + tris[i].p3[2]);
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_area(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed) lines = surf->mylines;
    else lines = surf->lines;
    double p12[3];
    for (int i = 0; i < nslocal; i++) {
      if (lines[i].mask & groupbit) {
	MathExtra::sub3(lines[i].p2,lines[i].p1,p12);
	buf[n] = MathExtra::len3(p12);
      } else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed) tris = surf->mytris;
    else tris = surf->tris;
    double p12[3],p23[3],cross[3];
    for (int i = 0; i < nslocal; i++) {
      if (tris[i].mask & groupbit) {
	MathExtra::sub3(tris[i].p2,tris[i].p1,p12);
	MathExtra::sub3(tris[i].p3,tris[i].p2,p12);
	MathExtra::cross3(p12,p23,cross);
	buf[n] = 0.5 * MathExtra::len3(cross);
      } else buf[n] = 0.0;
      n += nvalues;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_normx(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nslocal; i++) {
      if (lines[i].mask & groupbit) buf[n] = lines[i].norm[0];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nslocal; i++) {
      if (tris[i].mask & groupbit) buf[n] = tris[i].norm[0];
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_normy(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nslocal; i++) {
      if (lines[i].mask & groupbit) buf[n] = lines[i].norm[1];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nslocal; i++) {
      if (tris[i].mask & groupbit) buf[n] = tris[i].norm[1];
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurf::pack_normz(int n)
{
  Surf::Tri *tris;
  if (distributed) tris = surf->mytris;
  else tris = surf->tris;
  for (int i = 0; i < nslocal; i++) {
    if (tris[i].mask & groupbit) buf[n] = tris[i].norm[2];
    else buf[n] = 0.0;
    n += nvalues;
  }
}
