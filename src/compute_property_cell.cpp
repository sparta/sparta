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

#include "string.h"
#include "compute_property_cell.h"
#include "particle.h"
#include "grid.h"
#include "update.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

/* ---------------------------------------------------------------------- */

ComputePropertyCell::ComputePropertyCell(DSMC *dsmc, int narg, char **arg) :
  Compute(dsmc, narg, arg)
{
  if (narg < 2) error->all(FLERR,"Illegal compute property/cell command");

  per_cell_flag = 1;
  nvalues = narg - 2;
  if (nvalues == 1) size_per_cell_cols = 0;
  else size_per_cell_cols = nvalues;

  // parse input values
  // customize a new keyword by adding to if statement

  pack_choice = new FnPtrPack[nvalues];

  int i;
  int iarg = 2;
  while (iarg < narg) {
    i = iarg-2;

    if (strcmp(arg[iarg],"u") == 0) {
      pack_choice[i] = &ComputePropertyCell::pack_u;
    } else if (strcmp(arg[iarg],"v") == 0) {
      pack_choice[i] = &ComputePropertyCell::pack_v;
    } else if (strcmp(arg[iarg],"w") == 0) {
      pack_choice[i] = &ComputePropertyCell::pack_w;

    } else if (strcmp(arg[iarg],"usq") == 0) {
      pack_choice[i] = &ComputePropertyCell::pack_usq;
    } else if (strcmp(arg[iarg],"vsq") == 0) {
      pack_choice[i] = &ComputePropertyCell::pack_vsq;
    } else if (strcmp(arg[iarg],"wsq") == 0) {
      pack_choice[i] = &ComputePropertyCell::pack_wsq;

    } else error->all(FLERR,"Invalid keyword in compute property/cell command");
  }

  nglocal = grid->nlocal;
  vector = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputePropertyCell::~ComputePropertyCell()
{
  delete [] pack_choice;
}

/* ---------------------------------------------------------------------- */

void ComputePropertyCell::init()
{
}

/* ---------------------------------------------------------------------- */

void ComputePropertyCell::compute_per_cell()
{
  invoked_per_cell = update->ntimestep;

  // fill vector or array with per-atom values

  if (nvalues == 0) {
    pack_standard();
  } else if (nvalues == 1) {
    (this->*pack_choice[0])(0);
  } else {
    for (int n = 0; n < nvalues; n++)
      (this->*pack_choice[n])(n);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local cell-based array
------------------------------------------------------------------------- */

double ComputePropertyCell::memory_usage()
{
  double bytes;
  if (nvalues == 0) bytes = nglocal * sizeof(double);
  else bytes = nglocal*nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   pack standard quantities into local cell-based array
------------------------------------------------------------------------- */

void ComputePropertyCell::pack_standard()
{
  int icell,ilocal;
  double *v;

  Grid::OneCell *cells = grid->cells;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  for (int i = 0; i < nlocal; i++) {
    icell = particles[i].icell;
    ilocal = cells[icell].local;
    v = particles[i].v;

    array[ilocal][0] += 1.0;
    array[ilocal][1] += v[0];
    array[ilocal][2] += v[1];
    array[ilocal][3] += v[2];
    array[ilocal][4] += v[0]*v[0];
    array[ilocal][5] += v[1]*v[1];
    array[ilocal][6] += v[2]*v[2];
  }
}

/* ----------------------------------------------------------------------
   count particles per cell
------------------------------------------------------------------------- */

void ComputePropertyCell::pack_count()
{
  int icell,ilocal;

  Grid::OneCell *cells = grid->cells;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  for (int i = 0; i < nlocal; i++) {
    icell = particles[i].icell;
    ilocal = cells[icell].local;
    array[ilocal][0] += 1.0;
  }
}

/* ----------------------------------------------------------------------
   one method for every keyword compute property/cell can output
   the cell property is packed into buf starting at n with stride nvalues
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

void ComputePropertyCell::pack_u(int n)
{
  int icell,ilocal;
  double *v;

  Grid::OneCell *cells = grid->cells;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  for (int i = 0; i < nlocal; i++) {
    icell = particles[i].icell;
    ilocal = cells[icell].local;
    v = particles[i].v;
    array[ilocal][n] += v[0];
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyCell::pack_v(int n)
{
  int icell,ilocal;
  double *v;

  Grid::OneCell *cells = grid->cells;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  for (int i = 0; i < nlocal; i++) {
    icell = particles[i].icell;
    ilocal = cells[icell].local;
    v = particles[i].v;
    array[ilocal][n] += v[1];
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyCell::pack_w(int n)
{
  int icell,ilocal;
  double *v;

  Grid::OneCell *cells = grid->cells;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  for (int i = 0; i < nlocal; i++) {
    icell = particles[i].icell;
    ilocal = cells[icell].local;
    v = particles[i].v;
    array[ilocal][n] += v[2];
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyCell::pack_usq(int n)
{
  int icell,ilocal;
  double *v;

  Grid::OneCell *cells = grid->cells;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  for (int i = 0; i < nlocal; i++) {
    icell = particles[i].icell;
    ilocal = cells[icell].local;
    v = particles[i].v;
    array[ilocal][n] += v[0]*v[0];
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyCell::pack_vsq(int n)
{
  int icell,ilocal;
  double *v;

  Grid::OneCell *cells = grid->cells;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  for (int i = 0; i < nlocal; i++) {
    icell = particles[i].icell;
    ilocal = cells[icell].local;
    v = particles[i].v;
    array[ilocal][n] += v[1]*v[1];
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyCell::pack_wsq(int n)
{
  int icell,ilocal;
  double *v;

  Grid::OneCell *cells = grid->cells;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  for (int i = 0; i < nlocal; i++) {
    icell = particles[i].icell;
    ilocal = cells[icell].local;
    v = particles[i].v;
    array[ilocal][n] += v[2]*v[2];
  }
}
