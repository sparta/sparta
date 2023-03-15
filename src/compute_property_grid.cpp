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
#include "compute_property_grid.h"
#include "grid.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ComputePropertyGrid::ComputePropertyGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute property/grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  nvalues = narg - 3;

  // parse input values
  // customize a new keyword by adding to if statement

  pack_choice = new FnPtrPack[nvalues];
  index = new int[nvalues];

  int i;
  for (int iarg = 3; iarg < narg; iarg++) {
    i = iarg-3;

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[i] = &ComputePropertyGrid::pack_id;
      index[i] = 0;
    } else if (strcmp(arg[iarg],"proc") == 0) {
      pack_choice[i] = &ComputePropertyGrid::pack_proc;
      index[i] = 1;
    } else if (strcmp(arg[iarg],"xlo") == 0) {
      pack_choice[i] = &ComputePropertyGrid::pack_xlo;
      index[i] = 2;
    } else if (strcmp(arg[iarg],"ylo") == 0) {
      pack_choice[i] = &ComputePropertyGrid::pack_ylo;
      index[i] = 3;
    } else if (strcmp(arg[iarg],"zlo") == 0) {
      if (domain->dimension == 2)
        error->all(FLERR,
                   "Invalid compute property/grid field for 2d simulation");
      pack_choice[i] = &ComputePropertyGrid::pack_zlo;
      index[i] = 4;
    } else if (strcmp(arg[iarg],"xhi") == 0) {
      pack_choice[i] = &ComputePropertyGrid::pack_xhi;
      index[i] = 5;
    } else if (strcmp(arg[iarg],"yhi") == 0) {
      pack_choice[i] = &ComputePropertyGrid::pack_yhi;
      index[i] = 6;
    } else if (strcmp(arg[iarg],"zhi") == 0) {
      if (domain->dimension == 2)
        error->all(FLERR,
                   "Invalid compute property/grid field for 2d simulation");
      pack_choice[i] = &ComputePropertyGrid::pack_zhi;
      index[i] = 7;
    } else if (strcmp(arg[iarg],"xc") == 0) {
      pack_choice[i] = &ComputePropertyGrid::pack_xc;
      index[i] = 8;
    } else if (strcmp(arg[iarg],"yc") == 0) {
      pack_choice[i] = &ComputePropertyGrid::pack_yc;
      index[i] = 9;
    } else if (strcmp(arg[iarg],"zc") == 0) {
      if (domain->dimension == 2)
        error->all(FLERR,
                   "Invalid compute property/grid field for 2d simulation");
      pack_choice[i] = &ComputePropertyGrid::pack_zc;
      index[i] = 10;
    } else if (strcmp(arg[iarg],"vol") == 0) {
      pack_choice[i] = &ComputePropertyGrid::pack_vol;
      index[i] = 11;
    } else error->all(FLERR,"Invalid keyword in compute property/grid command");
  }

  per_grid_flag = 1;
  if (nvalues == 1) size_per_grid_cols = 0;
  else size_per_grid_cols = nvalues;

  nglocal = 0;
  vector_grid = NULL;
  array_grid = NULL;
}

/* ---------------------------------------------------------------------- */

ComputePropertyGrid::~ComputePropertyGrid()
{
  if (copymode) return;
  delete [] pack_choice;
  delete [] index;
  memory->destroy(vector_grid);
  memory->destroy(array_grid);
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::init()
{
  reallocate();
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  // fill vector or array with per-grid values

  if (nvalues == 1) {
    buf = vector_grid;
    (this->*pack_choice[0])(0);
  } else {
    if (nglocal) buf = &array_grid[0][0];
    else buf = NULL;
    for (int n = 0; n < nvalues; n++)
      (this->*pack_choice[n])(n);
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void ComputePropertyGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  nglocal = grid->nlocal;
  if (nvalues == 1) {
    memory->destroy(vector_grid);
    memory->create(vector_grid,nglocal,"property/grid:vector_grid");
  } else {
    memory->destroy(array_grid);
    memory->create(array_grid,nglocal,nvalues,"property/grid:array_grid");
  }
}

/* ----------------------------------------------------------------------
   memory usage of local grid-based array
------------------------------------------------------------------------- */

bigint ComputePropertyGrid::memory_usage()
{
  bigint bytes;
  bytes = nvalues*nglocal * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword compute property/grid can output
   the grid property is packed into buf starting at n with stride nvalues
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_id(int n)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  // NOTE: cellint (bigint) won't fit in double in some cases

  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].mask & groupbit) buf[n] = cells[i].id;
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_proc(int n)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].mask & groupbit) buf[n] = cells[i].proc;
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_xlo(int n)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].mask & groupbit) buf[n] = cells[i].lo[0];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_ylo(int n)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].mask & groupbit) buf[n] = cells[i].lo[1];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_zlo(int n)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].mask & groupbit) buf[n] = cells[i].lo[2];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_xhi(int n)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].mask & groupbit) buf[n] = cells[i].hi[0];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_yhi(int n)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].mask & groupbit) buf[n] = cells[i].hi[1];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_zhi(int n)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].mask & groupbit) buf[n] = cells[i].hi[2];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_xc(int n)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].mask & groupbit)
      buf[n] = 0.5 * (cells[i].lo[0] + cells[i].hi[0]);
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_yc(int n)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].mask & groupbit)
      buf[n] = 0.5 * (cells[i].lo[1] + cells[i].hi[1]);
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_zc(int n)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].mask & groupbit)
      buf[n] = 0.5 * (cells[i].lo[2] + cells[i].hi[2]);
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGrid::pack_vol(int n)
{
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].mask & groupbit) buf[n] = cinfo[i].volume;
    else buf[n] = 0.0;
    n += nvalues;
  }
}
