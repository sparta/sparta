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
#include "fix_temp_rescale.h"
#include "update.h"
#include "grid.h"
#include "particle.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixTempRescale::FixTempRescale(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix temp/rescale command");

  nevery = atoi(arg[2]);
  tstart = atof(arg[3]);
  tstop = atof(arg[4]);

  if (nevery <= 0) error->all(FLERR,"Illegal fix temp/rescale command");
  if (tstart < 0.0 || tstop < 0.0)
    error->all(FLERR,"Illegal fix temp/rescale command");

  // optional keyword

  aveflag = 0;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid fix temp/rescale command");
      if (strcmp(arg[iarg+1],"yes") == 0) aveflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) aveflag = 0;
      else error->all(FLERR,"Invalid fix temp/rescale command");
      iarg += 2;
    } else error->all(FLERR,"Invalid fix temp/rescale command");
  }

  // per-cell array for aveflag = 1 case

  maxgrid = 0;
  vcom = NULL;
}

/* ---------------------------------------------------------------------- */

FixTempRescale::~FixTempRescale()
{
  if (copymode) return;

  memory->destroy(vcom);
}

/* ---------------------------------------------------------------------- */

int FixTempRescale::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempRescale::init()
{
  tprefactor = update->mvv2e / (3.0*update->boltz);
}

/* ---------------------------------------------------------------------- */

void FixTempRescale::end_of_step()
{
  if (update->ntimestep % nevery) return;

  // set current t_target

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  double t_target = tstart + delta * (tstop-tstart);

  // sort particles by grid cell if needed

  if (!particle->sorted) particle->sort();

  // 2 variants of thermostatting

  if (!aveflag) end_of_step_no_average(t_target);
  else end_of_step_average(t_target);
}

/* ----------------------------------------------------------------------
   current thermal temperature is calculated on a per-cell basis
---------------------------------------------------------------------- */

void FixTempRescale::end_of_step_no_average(double t_target)
{
  // loop over grid cells and twice over particles in each cell
  // 1st pass: calc thermal temp via same logic as in ComputeThermalGrid
  // 2nd pass: rescale thermal velocity components

  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;
  int *next = particle->next;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  // loop over grid cells with more than 1 particle

  int ip,ispecies,count;
  double mass;
  double totmass,mvx,mvy,mvz,mvsq;
  double invtotmass,vscale;
  double vxcom,vycom,vzcom;
  double t_current;
  double *v;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cinfo[icell].count <= 1) continue;

    count = cinfo[icell].count;
    totmass = 0.0;
    mvx = mvy = mvz = 0.0;
    mvsq = 0.0;

    // 1st pass: loop over particles in cell
    // 5 tallies per particle: Mass, mVx, mVy, mVz, mV^2

    ip = cinfo[icell].first;
    while (ip >= 0) {
      ispecies = particles[ip].ispecies;
      mass = species[ispecies].mass;
      v = particles[ip].v;

      totmass += mass;
      mvx += mass*v[0];
      mvy += mass*v[1];
      mvz += mass*v[2];
      mvsq += mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

      ip = next[ip];
    }

    // COM velocity of particles in grid cell

    invtotmass = 1.0/totmass;
    vxcom = mvx * invtotmass;
    vycom = mvy * invtotmass;
    vzcom = mvz * invtotmass;

    // t_current = thermal T of particles in grid cell

    t_current = mvsq - (mvx*mvx + mvy*mvy + mvz*mvz)*invtotmass;
    t_current *= tprefactor/count;

    // vscale = scale factor for thermal velocity components

    vscale = sqrt(t_target/t_current);

    // 2nd pass: loop over particles in cell
    // rescale thermal velocity components

    ip = cinfo[icell].first;
    while (ip >= 0) {
      ispecies = particles[ip].ispecies;
      mass = species[ispecies].mass;
      v = particles[ip].v;

      v[0] = vscale*(v[0]-vxcom) + vxcom;
      v[1] = vscale*(v[1]-vycom) + vycom;
      v[2] = vscale*(v[2]-vzcom) + vzcom;

      ip = next[ip];
    }
  }
}

/* ----------------------------------------------------------------------
   current thermal temperature is averaged over all per-cell temperatures
---------------------------------------------------------------------- */

void FixTempRescale::end_of_step_average(double t_target)
{
  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;
  int *next = particle->next;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  // resize vcom if needed

  if (nglocal > maxgrid) {
    memory->destroy(vcom);
    maxgrid = nglocal + grid->nghost;
    memory->create(vcom,maxgrid,3,"temp/rescale:vcom");
  }

  // loop over grid cells to compute thermal T of each
  // only unsplit and sub cells, skip split cells

  int ip,ispecies,count;
  double mass;
  double totmass,mvx,mvy,mvz,mvsq;
  double invtotmass,t_one;
  double *v;

  bigint n_current_mine = 0.0;
  double t_current_mine = 0.0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit > 1) continue;

    count = cinfo[icell].count;
    totmass = 0.0;
    mvx = mvy = mvz = 0.0;
    mvsq = 0.0;

    // loop over particles in cell
    // 5 tallies per particle: Mass, mVx, mVy, mVz, mV^2

    ip = cinfo[icell].first;
    while (ip >= 0) {
      ispecies = particles[ip].ispecies;
      mass = species[ispecies].mass;
      v = particles[ip].v;

      totmass += mass;
      mvx += mass*v[0];
      mvy += mass*v[1];
      mvz += mass*v[2];
      mvsq += mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

      ip = next[ip];
    }

    // t_one = thermal T of particles in the cell
    // vcom = COM velocity of particle in the cell
    // if cell has <= 1 particle: set t_one = t_target, vcom = zero
    // likewise if t_one = 0.0: set t_one = t_target, vcom = zero
    //   corner case when all particles have same velocity

    if (count > 1) {
      invtotmass = 1.0/totmass;
      t_one = mvsq - (mvx*mvx + mvy*mvy + mvz*mvz)*invtotmass;
      t_one *= tprefactor/count;
      vcom[icell][0] = mvx * invtotmass;
      vcom[icell][1] = mvy * invtotmass;
      vcom[icell][2] = mvz * invtotmass;

    } else {
      t_one = t_target;
      vcom[icell][0] = vcom[icell][1] = vcom[icell][2] = 0.0;
    }

    if (t_one == 0.0) {
      t_one = t_target;
      vcom[icell][0] = vcom[icell][1] = vcom[icell][2] = 0.0;
    }

    // accumulate thermal T over my cells

    t_current_mine += t_one;
    n_current_mine++;
  }

  // t_current = average of cellwise thermal T across all cells
  // n_current = total # of cells contributing to t_current

  double t_current;
  MPI_Allreduce(&t_current_mine,&t_current,1,MPI_DOUBLE,MPI_SUM,world);

  bigint n_current;
  MPI_Allreduce(&n_current_mine,&n_current,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  // t_current = cellwise averaged thermal T
  // scale all particles in all cells by vscale

  t_current /= n_current;
  double vscale = sqrt(t_target/t_current);

  // loop over grid cells to rescale velocity of particles in each
  // single-particle cells are also rescaled, their vcom = 0.0

  for (int icell = 0; icell < nglocal; icell++) {
    if (cinfo[icell].count == 0) continue;

    // loop over particles in cell
    // rescale thermal velocity components

    ip = cinfo[icell].first;
    while (ip >= 0) {
      ispecies = particles[ip].ispecies;
      mass = species[ispecies].mass;
      v = particles[ip].v;

      v[0] = vscale*(v[0]-vcom[icell][0]) + vcom[icell][0];
      v[1] = vscale*(v[1]-vcom[icell][1]) + vcom[icell][1];
      v[2] = vscale*(v[2]-vcom[icell][2]) + vcom[icell][2];

      ip = next[ip];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double FixTempRescale::memory_usage()
{
  double bytes = 0.0;
  bytes += maxgrid*3 * sizeof(double);    // vcom
  return bytes;
}
