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
#include "create_molecules.h"
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

CreateMolecules::CreateMolecules(DSMC *dsmc) : Pointers(dsmc) {}

/* ---------------------------------------------------------------------- */

void CreateMolecules::command(int narg, char **arg)
{
  if (!domain->box_exist) 
    error->all(FLERR,
	       "Cannot create molecules before simulation box is defined");
  if (narg < 1) error->all(FLERR,"Illegal create_molecules command");

  imix = particle->find_mixture(arg[0]);
  if (imix < 0) error->all(FLERR,"Create_molecules mixture ID does not exist");
  particle->mixture[imix]->init();

  // optional args

  bigint np = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"n") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_molecules command");
      np = ATOBIGINT(arg[iarg+1]);
      if (np <= 0) error->all(FLERR,"Illegal create_molecules command");
      iarg += 2;
    } else error->all(FLERR,"Illegal create_molecules command");
  }

  // calculate Np if not set explicitly
  // NOTE: eventually adjust for cells with cut volume

  if (np == 0) {
    double voltotal;
    if (domain->dimension == 3)
      voltotal = domain->xprd * domain->yprd * domain->zprd;
    else voltotal = domain->xprd * domain->yprd;
    np = particle->mixture[imix]->nrho * voltotal / update->fnum;
  }

  // generate molecules

  bigint nprevious = particle->nglobal;
  create_local(np);

  // error check

  bigint nglobal;
  bigint nme = particle->nlocal;
  MPI_Allreduce(&nme,&nglobal,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  if (nglobal - nprevious != np) {
    char str[128];
    sprintf(str,"Created incorrect # of molecules: " 
	    BIGINT_FORMAT " versus " BIGINT_FORMAT,
	    nglobal-nprevious,np);
    error->all(FLERR,str);
  }
  particle->nglobal = nglobal;

  // print stats

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Created " BIGINT_FORMAT " molecules\n",np);
    if (logfile) fprintf(logfile,"Created " BIGINT_FORMAT " molecules\n",np);
  }
}

/* ----------------------------------------------------------------------
   create Np molecules in parallel
   every proc fraction of Np for cells it owns
   created particle attributes depend on number of procs
------------------------------------------------------------------------- */

void CreateMolecules::create_local(bigint np)
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
  // Nme = # of molecules I will create
  // MPI_Scan() logic insures sum of nme = Np
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
  memory->create(vols,nprocs,"create_molecules:vols");
  MPI_Allgather(&volupto,1,MPI_DOUBLE,vols,1,MPI_DOUBLE,world);

  bigint nstart,nstop;
  if (me > 0) nstart = np * vols[me-1]/vols[nprocs-1];
  else nstart = 0;
  nstop = np * vols[me]/vols[nprocs-1];
  bigint nme = nstop-nstart;

  memory->destroy(vols);

  // loop over cells I own
  // ntarget = floating point # of molecules to create in one cell
  // npercell = integer # of molecules to create in one cell
  // basing ntarget on accumulated volume and nprev insures Nme total creations
  // particle species = random value based on mixture fractions
  // particle velocity = stream velocity + thermal velocity

  int nspecies = particle->mixture[imix]->nspecies;
  int *species = particle->mixture[imix]->species;
  double *cummulative = particle->mixture[imix]->cummulative;
  double *vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;

  int ilocal,icell,npercell,ispecies;
  double x[3],v[3];
  double vol,ntarget,rn,vn,vr,theta1,theta2;
  double erote, ivib;

  double volsum = 0.0;
  bigint nprev = 0;

  for (int i = 0; i < nglocal; i++) {
    icell = mycells[i];
    lo = cells[icell].lo;
    hi = cells[icell].hi;

    if (dimension == 3) vol = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
    else vol = (hi[0]-lo[0]) * (hi[1]-lo[1]);
    volsum += vol;

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

      vn = vscale[ispecies] * sqrt(-log(random->uniform()));
      vr = vscale[ispecies] * sqrt(-log(random->uniform()));
      theta1 = MY_2PI * random->uniform();
      theta2 = MY_2PI * random->uniform();
/*    printf("%e %e %e %e\n", MY_2PI, vn, vr, vscale[ispecies]); */
	
      
      v[0] = vstream[0] + vn*cos(theta1);
      v[1] = vstream[1] + vr*sin(theta2);
      v[2] = vstream[2] + vr*sin(theta2);
/*
      erote = CreateMolecules.erot(isp);
      ivib = CreateMolecules.evib(isp);
*/
      particle->add_particle(0,ispecies,icell,x,v,erote,ivib);
      

    }

    nprev += npercell;
  }

  delete random;
}

/* ----------------------------------------------------------------------
   create Np molecules in serial
   every proc generates all Np coords, only keeps those in cells it owns
   created particle attributes should be independent of number of procs
------------------------------------------------------------------------- */

/*
void CreateMolecules::create_all(bigint n)
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

  // loop over all N molecules

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

double CreateMolecules::erot(int isp)
{
 RanPark *random = new RanPark(update->ranmaster->uniform());
 double erote,a,i,erm,b;

 if (particle->species[isp].rotdof == 2) {
  erote = -log(random->uniform()) * update->boltz * update->temp_thermal;
 }
 else {
  a=0.5*particle->species[isp].rotdof-1.;
  i=0;
  while (i == 0) {
    erm=random->uniform()*10.;
//-there is an energy cut-off at 10 kT
    b=pow(erm/a,a)*exp(a-erm);
    if (b > random->uniform()) i=1;
  }
    erote=erm*update->boltz*update->temp_thermal;
 }
  return erote;

}


int CreateMolecules::evib(int isp)
{
 RanPark *random = new RanPark(update->ranmaster->uniform());

 int ivib = -log(random->uniform()) * update->temp_thermal 
          / particle->species[isp].vibtemp;
 return ivib;

}

