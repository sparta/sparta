/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "fix_temp_rescale.h"
#include "update.h"
#include "grid.h"
#include "particle.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixTempRescale::FixTempRescale(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix temp/rescale command");

  nevery = atoi(arg[2]);
  tstart = atof(arg[3]);
  tstop = atof(arg[4]);

  if (nevery <= 0) error->all(FLERR,"Illegal fix temp/rescale command");
  if (tstart < 0.0 || tstop < 0.0)
    error->all(FLERR,"Illegal fix temp/rescale command");
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

  // loop over grid cells and twice over particles in each cell
  // 1st pass: calc thermal temp via same logic as in ComputeThermalGrid
  // 2nd pass: rescale thermal velocity components

  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;
  int *next = particle->next;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  // loop over grid cells with more than 1 particle

  int ip,ispecies;
  double mass;
  double count,totmass,mvx,mvy,mvz,mvsq;
  double invtotmass,vscale;
  double vxcom,vycom,vzcom;
  double t_current;
  double *v;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cinfo[icell].count <= 1) continue;

    count = 0.0;
    totmass = 0.0;
    mvx = mvy = mvz = 0.0;
    mvsq = 0.0;

    // 1st pass: loop over particles in cell
    // 6 tallies per particle: N, Mass, mVx, mVy, mVz, mV^2

    ip = cinfo[icell].first;
    while (ip >= 0) {
      ispecies = particles[ip].ispecies;
      mass = species[ispecies].mass;
      v = particles[ip].v;

      count += 1.0;
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
