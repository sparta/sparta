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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "create_particles.h"
#include "particle.h"
#include "update.h"
#include "grid.h"
#include "comm.h"
#include "domain.h"
#include "mixture.h"
#include "random_mars.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

CreateParticles::CreateParticles(DSMC *dsmc) : Pointers(dsmc) {}

/* ---------------------------------------------------------------------- */

void CreateParticles::command(int narg, char **arg)
{
  if (!domain->box_exist) 
    error->all(FLERR,
	       "Cannot create particles before simulation box is defined");
  if (narg < 1) error->all(FLERR,"Illegal create_particles command");

  imix = particle->find_mixture(arg[0]);
  if (imix < 0) error->all(FLERR,"Create_particles mixture ID does not exist");
  particle->mixture[imix]->init();

  // optional args

  bigint n = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"n") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_particles command");
      n = ATOBIGINT(arg[iarg+1]);
      if (n <= 0) error->all(FLERR,"Illegal create_particles command");
      iarg += 2;
    } else error->all(FLERR,"Illegal create_particles command");
  }

  // calculate N if not set explicitly
  // NOTE: eventually adjust for cells with cut volume

  if (n == 0) {
    double voltotal;
    if (domain->dimension == 3)
      voltotal = domain->xprd * domain->yprd * domain->zprd;
    else voltotal = domain->xprd * domain->yprd;
    n = update->nrho * voltotal / update->fnum;
  }

  // generate particles

  bigint nprevious = particle->nglobal;
  create_local(n);

  // error check

  bigint nglobal;
  bigint nme = particle->nlocal;
  MPI_Allreduce(&nme,&nglobal,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  if (nglobal - nprevious != n) {
    char str[128];
    sprintf(str,"Created incorrect # of particles: " 
	    BIGINT_FORMAT " versus " BIGINT_FORMAT,
	    nglobal-nprevious,n);
    error->all(FLERR,str);
  }
  particle->nglobal = nglobal;

  // print stats

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Created " BIGINT_FORMAT " particles\n",n);
    if (logfile) fprintf(logfile,"Created " BIGINT_FORMAT " particles\n",n);
  }
}

/* ----------------------------------------------------------------------
   create N particles in parallel
   every proc generates all N coords, only keeps those in cells it owns
   created particle attributes depend on number of procs
------------------------------------------------------------------------- */

void CreateParticles::create_local(bigint n)
{
  int dimension = domain->dimension;

  int me = comm->me;
  RanPark *random = new RanPark(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,me,100);

  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;
  int nglocal = grid->nlocal;

  // volme = volume of grid cells I own
  // Nme = # of particles I will create
  // MPI_Scan() logic insures sum of nme = N
  // NOTE: eventually adjust for cells with cut volume

  double *lo,*hi;
  double volme = 0.0;
  for (int i = 0; i < nglocal; i++) {
    lo = cells[mycells[i]].lo;
    hi = cells[mycells[i]].hi;
    if (dimension == 3) volme += (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
    else volme += (hi[0]-lo[0]) * (hi[1]-lo[1]);
  }
  
  double volupto;
  MPI_Scan(&volme,&volupto,1,MPI_DOUBLE,MPI_SUM,world);

  double *vols;
  int nprocs = comm->nprocs;
  memory->create(vols,nprocs,"create_particles:vols");
  MPI_Allgather(&volupto,1,MPI_DOUBLE,vols,1,MPI_DOUBLE,world);

  bigint nstart,nstop;
  if (me > 0) nstart = n * vols[me-1]/vols[nprocs-1];
  else nstart = 0;
  nstop = n * vols[me]/vols[nprocs-1];
  bigint nme = nstop-nstart;

  memory->destroy(vols);

  // loop over cells I own
  // ntarget = floating point # of particles to create in one cell
  // npercell = integer # of particles to create in one cell
  // basing ntarget on accumulated volume and nprev insures Nme total creations
  // particle species = random value based on mixture fractions
  // particle velocity = stream velocity + thermal velocity

  int nspecies = particle->mixture[imix]->nspecies;
  int *species = particle->mixture[imix]->species;
  double *cummulative = particle->mixture[imix]->cummulative;
  double **vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;

  int ilocal,icell,npercell,ispecies;
  double x[3],v[3];
  double volcell,ntarget,rn,vn,vr,theta1,theta2;

  double volsum = 0.0;
  bigint nprev = 0;

  for (int i = 0; i < nglocal; i++) {
    icell = mycells[i];
    lo = cells[icell].lo;
    hi = cells[icell].hi;

    if (dimension == 3) volcell = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
    else volcell = (hi[0]-lo[0]) * (hi[1]-lo[1]);
    volsum += volcell;

    ntarget = nme * volsum/volme - nprev;
    npercell = static_cast<int> (ntarget);
    if (random->uniform() < ntarget-npercell) npercell++;

    for (int m = 0; m < npercell; m++) {
      rn = random->uniform();
      ispecies = 0;
      while (cummulative[ispecies] < rn) ispecies++;

      x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
      x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
      x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
      if (dimension == 2) x[2] = 0.0;

      vn = vscale[ispecies] * random->gaussian();
      vr = vscale[ispecies] * random->gaussian();
      theta1 = MY_2PI * random->uniform();
      theta2 = MY_2PI * random->uniform();
	
      v[0] = vstream[ispecies][0] + vn*cos(theta1);
      v[1] = vstream[ispecies][1] + vr*sin(theta2);
      v[2] = vstream[ispecies][2] + vr*cos(theta2);

      particle->add_particle(0,ispecies,icell,x,v);
    }

    nprev += npercell;
  }

  delete random;
}

/* ----------------------------------------------------------------------
   create N particles in serial
   every proc generates all N coords, only keeps those in cells it owns
   created particle attributes should be independent of number of procs
------------------------------------------------------------------------- */

/*
void CreateParticles::create_all(bigint n)
{
  int dimension = domain->dimension;
  double xlo = domain->boxlo[0];
  double ylo = domain->boxlo[1];
  double zlo = domain->boxlo[2];
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  int me = comm->me;
  RanPark *random = new RandomPark(update->ranmaster->uniform());

  int icell;
  double x,y,z;

  // loop over all N particles

  for (bigint m = 0; m < n; m++) {
    x = xlo + random->uniform()*xprd;
    y = ylo + random->uniform()*yprd;
    z = zlo + random->uniform()*zprd;
    if (dimension == 2) z = 0.0;

    // which_cell() returns global grid cell index the particle is in
    // if I own that grid cell, add particle

    icell = grid->which_cell(x,y,z);
    if (grid->cells[icell].proc == me) {
      particle->add_particle(0,1,icell,x,y,z);
    }
  }

  delete random;
}
*/
