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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "create_particles.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "modify.h"
#include "comm.h"
#include "domain.h"
#include "region.h"
#include "random_mars.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};   // same as Grid

#define EPSZERO 1.0e-14

/* ---------------------------------------------------------------------- */

CreateParticles::CreateParticles(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

void CreateParticles::command(int narg, char **arg)
{
  if (!grid->exist)
    error->all(FLERR,"Cannot create particles before grid is defined");

  particle->exist = 1;

  if (narg < 1) error->all(FLERR,"Illegal create_particles command");

  imix = particle->find_mixture(arg[0]);
  if (imix < 0) error->all(FLERR,"Create_particles mixture ID does not exist");
  particle->mixture[imix]->init();

  // style arg

  bigint np = 0;
  single = 0;

  int iarg = 1;
  if (strcmp(arg[iarg],"n") == 0) {
    if (iarg+2 > narg) error->all(FLERR,"Illegal create_particles command");
    np = ATOBIGINT(arg[iarg+1]);
      if (np < 0) error->all(FLERR,"Illegal create_particles command");
      iarg += 2;
  } else if (strcmp(arg[iarg],"single") == 0) {
    if (iarg+8 > narg) error->all(FLERR,"Illegal create_particles command");
    single = 1;
    mspecies = particle->find_species(arg[iarg+1]);
    if (mspecies < 0) 
      error->all(FLERR,"Create_particles species ID does not exist");
    xp = atof(arg[iarg+2]);
    yp = atof(arg[iarg+3]);
    zp = atof(arg[iarg+4]);
    vx = atof(arg[iarg+5]);
    vy = atof(arg[iarg+6]);
    vz = atof(arg[iarg+7]);
    iarg += 8;
  } else error->all(FLERR,"Illegal create_particles command");

  // optional args

  int globalflag = 0;
  region = NULL;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"global") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_particles command");
      if (strcmp(arg[iarg+1],"no") == 0) globalflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) globalflag = 1;
      else error->all(FLERR,"Illegal create_particles command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_particles command");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion < 0) 
        error->all(FLERR,"Create_particles region does not exist");
      region = domain->regions[iregion];
      iarg += 2;
    } else error->all(FLERR,"Illegal create_particles command");
  }

  if (globalflag) 
    error->all(FLERR,"Create_particles global option not yet implemented");

  // calculate Np if not set explicitly

  if (single) np = 1;
  else if (np == 0) {
    Grid::ChildCell *cells = grid->cells;
    Grid::ChildInfo *cinfo = grid->cinfo;
    int nglocal = grid->nlocal;

    double flowvolme = 0.0;
    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit > 1) continue;
      if (cinfo[icell].type != INSIDE) 
        flowvolme += cinfo[icell].volume / cinfo[icell].weight;
    }
    double flowvol;
    MPI_Allreduce(&flowvolme,&flowvol,1,MPI_DOUBLE,MPI_SUM,world);
    np = particle->mixture[imix]->nrho * flowvol / update->fnum;
  }

  // generate particles
  // NOTE: invoke local or global option here

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  bigint nprevious = particle->nglobal;
  if (single) create_single();
  else if (!globalflag) create_local(np);
  //else create_global(np);

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // error check
  // only if no region specified

  bigint nglobal;
  bigint nme = particle->nlocal;
  MPI_Allreduce(&nme,&nglobal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  if (!region && nglobal-nprevious != np) {
    char str[128];
    sprintf(str,"Created incorrect # of particles: " 
	    BIGINT_FORMAT " versus " BIGINT_FORMAT,
	    nglobal-nprevious,np);
    error->all(FLERR,str);
  }
  bigint ncreated = nglobal-nprevious;
  particle->nglobal = nglobal;

  // print stats

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Created " BIGINT_FORMAT " particles\n",ncreated);
      fprintf(screen,"  CPU time = %g secs\n",time2-time1);
    }
    if (logfile) {
      fprintf(logfile,"Created " BIGINT_FORMAT " particles\n",ncreated);
      fprintf(logfile,"  CPU time = %g secs\n",time2-time1);
    }
  }
}

/* ----------------------------------------------------------------------
   create a single particle
   find cell it is in, and store on appropriate processor
------------------------------------------------------------------------- */

void CreateParticles::create_single()
{
  int i,m;
  double x[3],v[3],vstream[3];
  double *lo,*hi;

  x[0] = xp;  x[1] = yp;  x[2] = zp;
  v[0] = vx;  v[1] = vy;  v[2] = vz;
  double temp_thermal = particle->mixture[imix]->temp_thermal;
  vstream[0] = vstream[1] = vstream[2] = 0.0;

  if (domain->dimension == 2 && x[2] != 0.0)
    error->all(FLERR,"Create_particles single requires z = 0 "
	       "for 2d simulation");

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  int iwhich = -1;
  for (int icell = 0; icell < nglocal; icell++) {
    if (cinfo[icell].type != OUTSIDE) continue;
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    if (x[0] >= lo[0] && x[0] < hi[0] &&
	x[1] >= lo[1] && x[1] < hi[1] &&
	x[2] >= lo[2] && x[2] < hi[2]) iwhich = icell;
  }

  // insure that exactly one proc found cell to insert particle into

  int flag,flagall;
  if (iwhich < 0) flag = 0;
  else flag = 1;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall != 1) 
    error->all(FLERR,"Could not create a single particle");

  // nfix_add_particle = # of fixes with add_particle() method

  modify->list_init_fixes();
  int nfix_add_particle = modify->n_add_particle;

  // add the particle

  RanPark *random = new RanPark(update->ranmaster->uniform());

  if (iwhich >= 0) {
    int id = MAXSMALLINT*random->uniform();
    double erot = particle->erot(mspecies,temp_thermal,random);
    double evib = particle->evib(mspecies,temp_thermal,random);
    particle->add_particle(id,mspecies,iwhich,x,v,erot,evib);
    if (nfix_add_particle) 
      modify->add_particle(particle->nlocal-1,temp_thermal,vstream);
  }

  delete random;
}

/* ----------------------------------------------------------------------
   create Np particles in parallel
   every proc creates fraction of Np for cells it owns
   only insert in cells uncut by surfs
   account for cell weighting
   attributes of created particle depend on number of procs
------------------------------------------------------------------------- */

void CreateParticles::create_local(bigint np)
{
  int dimension = domain->dimension;

  int me = comm->me;
  RanPark *random = new RanPark(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,me,100);

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  // volme = volume of grid cells I own that are OUTSIDE surfs
  // skip cells entirely outside region
  // Nme = # of particles I will create
  // MPI_Scan() logic insures sum of nme = Np

  double *lo,*hi;
  double volone;

  double volme = 0.0;
  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].type != OUTSIDE) continue;
    lo = cells[i].lo;
    hi = cells[i].hi;
    if (region && region->bboxflag && outside_region(dimension,lo,hi))
      continue;

    if (dimension == 3) volone = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
    else if (domain->axisymmetric) 
      volone = (hi[0]-lo[0]) * (hi[1]*hi[1]-lo[1]*lo[1])*MY_PI;
    else volone = (hi[0]-lo[0]) * (hi[1]-lo[1]);
    volme += volone / cinfo[i].weight;
  }
  
  double volupto;
  MPI_Scan(&volme,&volupto,1,MPI_DOUBLE,MPI_SUM,world);

  double *vols;
  int nprocs = comm->nprocs;
  memory->create(vols,nprocs,"create_particles:vols");
  MPI_Allgather(&volupto,1,MPI_DOUBLE,vols,1,MPI_DOUBLE,world);

  // nme = # of particles for me to create
  // gathered Scan results not guaranteed to be monotonically increasing
  // can cause epsilon mis-counts for huge particle counts
  // enforce that by brute force

  for (int i = 1; i < nprocs; i++)
    if (vols[i] != vols[i-1] && 
        fabs(vols[i]-vols[i-1])/vols[nprocs-1] < EPSZERO)
      vols[i] = vols[i-1];

  bigint nstart,nstop;
  if (me > 0) nstart = static_cast<bigint> (np * (vols[me-1]/vols[nprocs-1]));
  else nstart = 0;
  nstop = static_cast<bigint> (np * (vols[me]/vols[nprocs-1]));
  bigint nme = nstop-nstart;

  memory->destroy(vols);

  // nfix_add_particle = # of fixes with add_particle() method

  modify->list_init_fixes();
  int nfix_add_particle = modify->n_add_particle;

  // loop over cells I own
  // only add particles to OUTSIDE cells
  // skip cells entirely outside region
  // ntarget = floating point # of particles to create in one cell
  // npercell = integer # of particles to create in one cell
  // basing ntarget on accumulated volume and nprev insures Nme total creations
  // particle species = random value based on mixture fractions
  // particle velocity = stream velocity + thermal velocity

  int nspecies = particle->mixture[imix]->nspecies;
  int *species = particle->mixture[imix]->species;
  double *cummulative = particle->mixture[imix]->cummulative;
  double *vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;
  double temp_thermal = particle->mixture[imix]->temp_thermal;

  int npercell,isp,ispecies,id;
  double x[3],v[3];
  double ntarget,rn,vn,vr,theta1,theta2,erot,evib;

  double volsum = 0.0;
  bigint nprev = 0;

  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].type != OUTSIDE) continue;
    lo = cells[i].lo;
    hi = cells[i].hi;
    if (region && region->bboxflag && outside_region(dimension,lo,hi))
      continue;

    if (dimension == 3) volone = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
    else if (domain->axisymmetric)
      volone = (hi[0]-lo[0]) * (hi[1]*hi[1]-lo[1]*lo[1])*MY_PI;
    else volone = (hi[0]-lo[0]) * (hi[1]-lo[1]);
    volsum += volone / cinfo[i].weight;

    ntarget = nme * volsum/volme - nprev;
    npercell = static_cast<int> (ntarget);
    if (random->uniform() < ntarget-npercell) npercell++;

    for (int m = 0; m < npercell; m++) {
      rn = random->uniform();
      isp = 0;
      while (cummulative[isp] < rn) isp++;
      ispecies = species[isp];

      x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
      x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
      x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
      if (dimension == 2) x[2] = 0.0;

      if (region && !region->match(x)) continue;

      // double boundary = 1.5E-3-1.E-4*sin((x[0]/0.5E-3)*2.*MY_PI+MY_PI*0.5);
      // double boundary = 1.5E-3-1.E-4*sin(x[2]/0.5E-3*2.*MY_PI+MY_PI*0.5)*
      //                   sin(x[0]/0.5E-3*2.*MY_PI+MY_PI*0.5);
      // if (x[1]>=boundary) ispecies = 1;
      // if (x[1]<boundary) ispecies = 0;

      vn = vscale[isp] * sqrt(-log(random->uniform()));
      vr = vscale[isp] * sqrt(-log(random->uniform()));
      theta1 = MY_2PI * random->uniform();
      theta2 = MY_2PI * random->uniform();
	
      v[0] = vstream[0] + vn*cos(theta1);
      v[1] = vstream[1] + vr*cos(theta2);
      v[2] = vstream[2] + vr*sin(theta2);

      erot = particle->erot(ispecies,temp_thermal,random);
      evib = particle->evib(ispecies,temp_thermal,random);

      id = MAXSMALLINT*random->uniform();

      particle->add_particle(id,ispecies,i,x,v,erot,evib);
      if (nfix_add_particle) 
        modify->add_particle(particle->nlocal-1,temp_thermal,vstream);
    }

    nprev += npercell;
  }

  delete random;
}

/* ----------------------------------------------------------------------
   return 1 if grid cell with lo/hi is entirely outside region bounding box
   else return 0
------------------------------------------------------------------------- */

int CreateParticles::outside_region(int dim, double *lo, double *hi)
{
  int flag = 1;
  if (hi[0] > region->extent_xlo &&
      lo[0] < region->extent_xhi) flag = 0;
  if (hi[1] > region->extent_ylo &&
      lo[1] < region->extent_yhi) flag = 0;
  if (dim == 3) {
    if (hi[2] > region->extent_zlo &&
        lo[2] < region->extent_zhi) flag = 0;
  }
  return flag;
}

/* ----------------------------------------------------------------------
   create Np particles in serial
   every proc generates all Np coords, only keeps those in cells it owns
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

  int icell,id;
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
    id = MAXSMALLINT*random->uniform();

    if (grid->cells[icell].proc == me) {
      particle->add_particle(id,1,icell,x,y,z);
    }
  }

  delete random;
}
*/
