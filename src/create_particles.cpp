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
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};   // same as Grid

#define MAXATTEMPT 1024      // max attempts to insert a particle into cut/split cell
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

  np = 0;
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

  cutflag = 1;
  int globalflag = 0;
  twopass = 0;
  region = NULL;
  speciesflag = densflag = velflag = tempflag = 0;
  sstr = sxstr = systr = szstr = NULL;
  dstr = dxstr = dystr = dzstr = NULL;
  tstr = txstr = tystr = tzstr = NULL;
  vxstr = vystr = vzstr = vstrx = vstry = vstrz = NULL;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"cut") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_particles command");
      if (strcmp(arg[iarg+1],"no") == 0) cutflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) cutflag = 1;
      else error->all(FLERR,"Illegal create_particles command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"global") == 0) {
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
    } else if (strcmp(arg[iarg],"species") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal create_particles command");
      speciesflag = 1;
      sstr = arg[iarg+1];
      if (strcmp(arg[iarg+2],"NULL") == 0) sxstr = NULL;
      else sxstr = arg[iarg+2];
      if (strcmp(arg[iarg+3],"NULL") == 0) systr = NULL;
      else systr = arg[iarg+3];
      if (strcmp(arg[iarg+4],"NULL") == 0) szstr = NULL;
      else szstr = arg[iarg+4];
      iarg += 5;
    } else if (strcmp(arg[iarg],"density") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal create_particles command");
      densflag = 1;
      dstr = arg[iarg+1];
      if (strcmp(arg[iarg+2],"NULL") == 0) dxstr = NULL;
      else dxstr = arg[iarg+2];
      if (strcmp(arg[iarg+3],"NULL") == 0) dystr = NULL;
      else dystr = arg[iarg+3];
      if (strcmp(arg[iarg+4],"NULL") == 0) dzstr = NULL;
      else dzstr = arg[iarg+4];
      iarg += 5;
    } else if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal create_particles command");
      tempflag = 1;
      tstr = arg[iarg+1];
      if (strcmp(arg[iarg+2],"NULL") == 0) txstr = NULL;
      else txstr = arg[iarg+2];
      if (strcmp(arg[iarg+3],"NULL") == 0) tystr = NULL;
      else tystr = arg[iarg+3];
      if (strcmp(arg[iarg+4],"NULL") == 0) tzstr = NULL;
      else tzstr = arg[iarg+4];
      iarg += 5;
    } else if (strcmp(arg[iarg],"velocity") == 0) {
      if (iarg+7 > narg) error->all(FLERR,"Illegal create_particles command");
      velflag = 1;
      if (strcmp(arg[iarg+1],"NULL") == 0) vxstr = NULL;
      else vxstr = arg[iarg+1];
      if (strcmp(arg[iarg+2],"NULL") == 0) vystr = NULL;
      else vystr = arg[iarg+2];
      if (strcmp(arg[iarg+3],"NULL") == 0) vzstr = NULL;
      else vzstr = arg[iarg+3];
      if (strcmp(arg[iarg+4],"NULL") == 0) vstrx = NULL;
      else vstrx = arg[iarg+4];
      if (strcmp(arg[iarg+5],"NULL") == 0) vstry = NULL;
      else vstry = arg[iarg+5];
      if (strcmp(arg[iarg+6],"NULL") == 0) vstrz = NULL;
      else vstrz = arg[iarg+6];
      iarg += 7;
    } else if (strcmp(arg[iarg],"twopass") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal create_particles command");
      twopass = 1;
      iarg += 1;
    } else error->all(FLERR,"Illegal create_particles command");
  }

  if (globalflag)
    error->all(FLERR,"Create_particles global option not yet implemented");

  // error checks and further setup for variables

  if (speciesflag) {
    svar = input->variable->find(sstr);
    if (svar < 0)
      error->all(FLERR,"Variable name for create_particles does not exist");
    if (!input->variable->equal_style(svar))
      error->all(FLERR,"Variable for create_particles is invalid style");
    if (sxstr) {
      sxvar = input->variable->find(sxstr);
      if (sxvar < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->internal_style(sxvar))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
    if (systr) {
      syvar = input->variable->find(systr);
      if (syvar < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->internal_style(syvar))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
    if (szstr) {
      szvar = input->variable->find(szstr);
      if (szvar < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->internal_style(szvar))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
  }

  if (densflag) {
    dvar = input->variable->find(dstr);
    if (dvar < 0)
      error->all(FLERR,"Variable name for create_particles does not exist");
    if (!input->variable->equal_style(dvar))
      error->all(FLERR,"Variable for create_particles is invalid style");
    if (dxstr) {
      dxvar = input->variable->find(dxstr);
      if (dxvar < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->internal_style(dxvar))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
    if (dystr) {
      dyvar = input->variable->find(dystr);
      if (dyvar < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->internal_style(dyvar))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
    if (dzstr) {
      dzvar = input->variable->find(dzstr);
      if (dzvar < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->internal_style(dzvar))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
  }

  if (tempflag) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR,"Variable name for create_particles does not exist");
    if (!input->variable->equal_style(tvar))
      error->all(FLERR,"Variable for create_particles is invalid style");
    if (txstr) {
      txvar = input->variable->find(txstr);
      if (txvar < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->internal_style(txvar))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
    if (tystr) {
      tyvar = input->variable->find(tystr);
      if (tyvar < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->internal_style(tyvar))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
    if (tzstr) {
      tzvar = input->variable->find(tzstr);
      if (tzvar < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->internal_style(tzvar))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
  }

  if (velflag) {
    if (vxstr) {
      vxvar = input->variable->find(vxstr);
      if (vxvar < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->equal_style(vxvar))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
    if (vystr) {
      vyvar = input->variable->find(vystr);
      if (vyvar < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->equal_style(vyvar))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
    if (vzstr) {
      vzvar = input->variable->find(vzstr);
      if (vzvar < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->equal_style(vzvar))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
    if (vstrx) {
      vvarx = input->variable->find(vstrx);
      if (vvarx < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->internal_style(vvarx))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
    if (vstry) {
      vvary = input->variable->find(vstry);
      if (vvary < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->internal_style(vvary))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
    if (vstrz) {
      vvarz = input->variable->find(vstrz);
      if (vvarz < 0)
        error->all(FLERR,"Variable name for create_particles does not exist");
      if (!input->variable->internal_style(vvarz))
        error->all(FLERR,"Variable for create_particles is invalid style");
    }
  }

  // generate particles

  if (comm->me == 0)
    if (screen) fprintf(screen,"Creating particles ...\n");

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  bigint nprevious = particle->nglobal;

  if (single) create_single();
  else if (!globalflag) {
    if (twopass) create_local_twopass();
    else create_local();
  } else {
    error->all(FLERR,"Create_particles global option not yet implemented");
    // create_global();
  }

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // issue warning if created particle count is unexpected
  // only if no region and no variable density specified

  bigint nglobal;
  bigint nme = particle->nlocal;
  MPI_Allreduce(&nme,&nglobal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  if (!region && !densflag && nglobal-nprevious != np) {
    char str[128];
    sprintf(str,"Created unexpected # of particles: "
	    BIGINT_FORMAT " versus " BIGINT_FORMAT,
	    nglobal-nprevious,np);
    if (comm->me == 0) error->warning(FLERR,str);
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
  np = 1;

  double x[3],v[3],vstream[3];
  double *lo,*hi;

  x[0] = xp;  x[1] = yp;  x[2] = zp;
  v[0] = vx;  v[1] = vy;  v[2] = vz;
  double temp_thermal = particle->mixture[imix]->temp_thermal;
  double temp_rot = particle->mixture[imix]->temp_rot;
  double temp_vib = particle->mixture[imix]->temp_vib;
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

  // nfix_update_custom = # of fixes with update_custom() method

  particle->error_custom();
  modify->list_init_fixes();
  int nfix_update_custom = modify->n_update_custom;

  // add the particle

  RanKnuth *random = new RanKnuth(update->ranmaster->uniform());

  if (iwhich >= 0) {
    int id = MAXSMALLINT*random->uniform();
    double erot = particle->erot(mspecies,temp_rot,random);
    double evib = particle->evib(mspecies,temp_vib,random);
    particle->add_particle(id,mspecies,iwhich,x,v,erot,evib);
    if (nfix_update_custom)
      modify->update_custom(particle->nlocal-1,temp_thermal,
                           temp_rot,temp_vib,vstream);
  }

  delete random;
}

/* ----------------------------------------------------------------------
   create particles in parallel
   every proc creates fraction of particles for cells it owns
   cutflag determines whether to insert in all cells or only ones uncut by surfs
   account for cell weighting
   attributes of created particle depend on number of procs
------------------------------------------------------------------------- */

void CreateParticles::create_local()
{
  int dimension = domain->dimension;

  int me = comm->me;
  RanKnuth *random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,me,100);

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::SplitInfo *sinfo = grid->sinfo;
  int nglocal = grid->nlocal;

  // flowvol = total weighted flow volume of all cells
  //   skip cells inside surfs and split cells
  //   skip cells outside defined region
  // insertvol = subset of flowvol for cells eligible for insertion
  //   insertvol = flowvol if cutflag = 1
  //   insertvol < flowvol possible if cutflag = 0 (no cut cells)

  double flowvolme = 0.0;
  double insertvolme = 0.0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cinfo[icell].type == INSIDE) continue;
    if (cells[icell].nsplit > 1) continue;
    if (region && region->bboxflag &&
        outside_region(dimension,cells[icell].lo,cells[icell].hi))
      continue;

    flowvolme += cinfo[icell].volume / cinfo[icell].weight;
    if (!cutflag && cells[icell].nsurf) continue;
    insertvolme += cinfo[icell].volume / cinfo[icell].weight;
  }

  // calculate total Np if not set explicitly
  // based on total flowvol and mixture density

  if (np == 0) {
    double flowvol;
    MPI_Allreduce(&flowvolme,&flowvol,1,MPI_DOUBLE,MPI_SUM,world);
    np = particle->mixture[imix]->nrho * flowvol / update->fnum;
  }

  // gather cummulative insertion volumes across all procs

  double volupto;
  MPI_Scan(&insertvolme,&volupto,1,MPI_DOUBLE,MPI_SUM,world);

  double *vols;
  int nprocs = comm->nprocs;
  memory->create(vols,nprocs,"create_particles:vols");
  MPI_Allgather(&volupto,1,MPI_DOUBLE,vols,1,MPI_DOUBLE,world);

  // gathered Scan results not guaranteed to be monotonically increasing
  //   can cause epsilon mis-counts for huge particle counts
  //   so enforce monotonic increase by brute force

  for (int i = 1; i < nprocs; i++)
    if (vols[i] != vols[i-1] &&
        fabs(vols[i]-vols[i-1])/vols[nprocs-1] < EPSZERO)
      vols[i] = vols[i-1];

  // nme = # of particles for me to create
  // based on fraction of insertvol I own
  // loop over procs insures sum of nme = Np

  bigint nstart,nstop;
  if (me > 0) nstart = static_cast<bigint> (np * (vols[me-1]/vols[nprocs-1]));
  else nstart = 0;
  nstop = static_cast<bigint> (np * (vols[me]/vols[nprocs-1]));
  bigint nme = nstop-nstart;

  memory->destroy(vols);

  // nfix_update_custom = # of fixes with update_custom() method

  particle->error_custom();
  modify->list_init_fixes();
  int nfix_update_custom = modify->n_update_custom;

  // loop over cells I own
  // only add particles to cells eligible for insertion
  // ntarget = floating point # of particles to create in one cell
  // npercell = integer # of particles to create in one cell
  // basing ntarget on accumulated volume and nprev insures Nme total creations
  // particle species = random value based on mixture fractions
  // particle velocity = stream velocity + thermal velocity

  int *species = particle->mixture[imix]->species;
  double *cummulative = particle->mixture[imix]->cummulative;
  double *vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;
  int nspecies = particle->mixture[imix]->nspecies;
  double temp_thermal = particle->mixture[imix]->temp_thermal;
  double temp_rot = particle->mixture[imix]->temp_rot;
  double temp_vib = particle->mixture[imix]->temp_vib;

  int npercell,ncreate,isp,ispecies,id,pflag,subcell;
  double x[3],v[3],xcell[3],vstream_variable[3];
  double ntarget,scale,rn,vn,vr,theta1,theta2,erot,evib;
  double *lo,*hi;

  double tempscale = 1.0;
  double sqrttempscale = 1.0;

  double volsum = 0.0;
  bigint nprev = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cinfo[icell].type == INSIDE) continue;
    if (cells[icell].nsplit > 1) continue;
    if (cinfo[icell].volume == 0.0) continue;
    if (region && region->bboxflag &&
        outside_region(dimension,cells[icell].lo,cells[icell].hi))
      continue;
    if (!cutflag && cells[icell].nsurf) continue;

    volsum += cinfo[icell].volume / cinfo[icell].weight;

    ntarget = nme * volsum/insertvolme - nprev;
    npercell = static_cast<int> (ntarget);

    if (random->uniform() < ntarget-npercell) npercell++;
    ncreate = npercell;

    lo = cells[icell].lo;
    hi = cells[icell].hi;

    if (densflag) {
      scale = density_variable(lo,hi);
      ntarget *= scale;
      ncreate = static_cast<int> (ntarget);
      if (random->uniform() < ntarget-ncreate) ncreate++;
    }

    // if surfs in cell, use xcell for all created particle attempts

    if (cells[icell].nsurf)
      pflag = grid->point_outside_surfs(icell,xcell);

    for (int m = 0; m < ncreate; m++) {

      // generate random position X for new particle

      x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
      x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
      x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
      if (dimension == 2) x[2] = 0.0;

      // if surfs, check if random position is in flow region
      // if subcell, also check if in correct subcell
      // if not, attempt new insertion up to MAXATTEMPT times

      if (cells[icell].nsurf && pflag) {
        int nattempt = 0;
        while (nattempt < MAXATTEMPT) {
          if (grid->outside_surfs(icell,x,xcell)) {
            if (cells[icell].nsplit == 1) break;
            int splitcell = sinfo[cells[icell].isplit].icell;
            if (dimension == 2) subcell = update->split2d(splitcell,x);
            else subcell = update->split3d(splitcell,x);
            if (subcell == icell) break;
          }

          nattempt++;

          x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
          x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
          x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
          if (dimension == 2) x[2] = 0.0;
        }

        // particle insertion was unsuccessful

        if (nattempt >= MAXATTEMPT) continue;
      }

      // if region defined, skip if particle outside region

      if (region && !region->match(x)) continue;

      // insertion of particle at position X is accepted
      // calculate all other particle properties

      rn = random->uniform();

      isp = 0;
      while (cummulative[isp] < rn) isp++;
      ispecies = species[isp];

      if (speciesflag) {
        isp = species_variable(x) - 1;
        if (isp < 0 || isp >= nspecies) continue;
        ispecies = species[isp];
      }

      if (tempflag) {
        tempscale = temperature_variable(x);
        sqrttempscale = sqrt(tempscale);
      }

      vn = vscale[isp] * sqrttempscale * sqrt(-log(random->uniform()));
      vr = vscale[isp] * sqrttempscale * sqrt(-log(random->uniform()));
      theta1 = MY_2PI * random->uniform();
      theta2 = MY_2PI * random->uniform();

      if (velflag) {
        velocity_variable(x,vstream,vstream_variable);
        v[0] = vstream_variable[0] + vn*cos(theta1);
        v[1] = vstream_variable[1] + vr*cos(theta2);
        v[2] = vstream_variable[2] + vr*sin(theta2);
      } else {
        v[0] = vstream[0] + vn*cos(theta1);
        v[1] = vstream[1] + vr*cos(theta2);
        v[2] = vstream[2] + vr*sin(theta2);
      }

      erot = particle->erot(ispecies,temp_rot*tempscale,random);
      evib = particle->evib(ispecies,temp_vib*tempscale,random);

      id = MAXSMALLINT*random->uniform();

      particle->add_particle(id,ispecies,icell,x,v,erot,evib);

      if (nfix_update_custom)
        modify->update_custom(particle->nlocal-1,temp_thermal,
                             temp_rot,temp_vib,vstream);
    }

    // increment count without effect of density variation
    // so that target insertion count is undisturbed

    nprev += npercell;
  }

  delete random;
}

/* ----------------------------------------------------------------------
   create particles in parallel in two passes
   every proc creates fraction of paricles for cells it owns
   cutflag determines whether to insert in all cells or only ones uncut by surfs
   account for cell weighting
   attributes of created particle depend on number of procs
   use this version to be compatible with how Kokkos creates particles
------------------------------------------------------------------------- */

void CreateParticles::create_local_twopass()
{
  int dimension = domain->dimension;

  int me = comm->me;
  RanKnuth *random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,me,100);

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::SplitInfo *sinfo = grid->sinfo;
  int nglocal = grid->nlocal;

  // flowvol = total weighted flow volume of all cells
  //   skip cells inside surfs and split cells
  //   skip cells outside defined region
  // insertvol = subset of flowvol for cells eligible for insertion
  //   insertvol = flowvol if cutflag = 1
  //   insertvol < flowvol possible if cutflag = 0 (no cut cells)

  double flowvolme = 0.0;
  double insertvolme = 0.0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cinfo[icell].type == INSIDE) continue;
    if (cells[icell].nsplit > 1) continue;
    if (region && region->bboxflag &&
        outside_region(dimension,cells[icell].lo,cells[icell].hi))
      continue;

    flowvolme += cinfo[icell].volume / cinfo[icell].weight;
    if (!cutflag && cells[icell].nsurf) continue;
    insertvolme += cinfo[icell].volume / cinfo[icell].weight;
  }

  // calculate total Np if not set explicitly
  // based on total flowvol and mixture density

  if (np == 0) {
    double flowvol;
    MPI_Allreduce(&flowvolme,&flowvol,1,MPI_DOUBLE,MPI_SUM,world);
    np = particle->mixture[imix]->nrho * flowvol / update->fnum;
  }

  // gather cummulative insertion volumes across all procs

  double volupto;
  MPI_Scan(&insertvolme,&volupto,1,MPI_DOUBLE,MPI_SUM,world);

  double *vols;
  int nprocs = comm->nprocs;
  memory->create(vols,nprocs,"create_particles:vols");
  MPI_Allgather(&volupto,1,MPI_DOUBLE,vols,1,MPI_DOUBLE,world);

  // gathered Scan results not guaranteed to be monotonically increasing
  //   can cause epsilon mis-counts for huge particle counts
  //   so enforce monotonic increase by brute force

  for (int i = 1; i < nprocs; i++)
    if (vols[i] != vols[i-1] &&
        fabs(vols[i]-vols[i-1])/vols[nprocs-1] < EPSZERO)
      vols[i] = vols[i-1];

  // nme = # of particles for me to create
  // based on fraction of insertvol I own
  // loop over procs insures sum of nme = Np

  bigint nstart,nstop;
  if (me > 0) nstart = static_cast<bigint> (np * (vols[me-1]/vols[nprocs-1]));
  else nstart = 0;
  nstop = static_cast<bigint> (np * (vols[me]/vols[nprocs-1]));
  bigint nme = nstop-nstart;

  memory->destroy(vols);

  // nfix_update_custom = # of fixes with update_custom() method

  modify->list_init_fixes();
  int nfix_update_custom = modify->n_update_custom;

  // loop over cells I own
  // only add particles to cells eligible for insertion
  // ntarget = floating point # of particles to create in one cell
  // npercell = integer # of particles to create in one cell
  // basing ntarget on accumulated volume and nprev insures Nme total creations
  // particle species = random value based on mixture fractions
  // particle velocity = stream velocity + thermal velocity

  int *species = particle->mixture[imix]->species;
  double *cummulative = particle->mixture[imix]->cummulative;
  double *vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;
  int nspecies = particle->mixture[imix]->nspecies;
  double temp_thermal = particle->mixture[imix]->temp_thermal;
  double temp_rot = particle->mixture[imix]->temp_rot;
  double temp_vib = particle->mixture[imix]->temp_vib;

  int npercell,ncreate,isp,ispecies,id,pflag,subcell;
  double x[3],v[3],xcell[3],vstream_variable[3];
  double ntarget,scale,rn,vn,vr,theta1,theta2,erot,evib;
  double *lo,*hi;

  double tempscale = 1.0;
  double sqrttempscale = 1.0;

  double volsum = 0.0;
  bigint nprev = 0;

  // first pass, just calculate # of particles to create
  // ncreate_values[icell] = # of particles to create in ICELL

  int *ncreate_values;
  memory->create(ncreate_values, nglocal, "create_particles:ncreate");

  for (int icell = 0; icell < nglocal; icell++) {
    if (cinfo[icell].type == INSIDE) continue;
    if (cells[icell].nsplit > 1) continue;
    if (cinfo[icell].volume == 0.0) continue;
    if (region && region->bboxflag &&
        outside_region(dimension,cells[icell].lo,cells[icell].hi))
      continue;
    if (!cutflag && cells[icell].nsurf) continue;

    volsum += cinfo[icell].volume / cinfo[icell].weight;

    ntarget = nme * volsum/insertvolme - nprev;
    npercell = static_cast<int> (ntarget);

    if (random->uniform() < ntarget-npercell) npercell++;
    ncreate = npercell;

    lo = cells[icell].lo;
    hi = cells[icell].hi;

    if (densflag) {
      scale = density_variable(lo,hi);
      ntarget *= scale;
      ncreate = static_cast<int> (ntarget);
      if (random->uniform() < ntarget-ncreate) ncreate++;
    }

    ncreate_values[icell] = ncreate;

    // increment count without effect of density variation
    // so that target insertion count is undisturbed

    nprev += npercell;
  }

  // second pass, create particles using ncreate_values

  for (int icell = 0; icell < nglocal; icell++) {
    if (cinfo[icell].type == INSIDE) continue;
    if (cells[icell].nsplit > 1) continue;
    if (cinfo[icell].volume == 0.0) continue;
    if (region && region->bboxflag &&
        outside_region(dimension,cells[icell].lo,cells[icell].hi))
      continue;
    if (!cutflag && cells[icell].nsurf) continue;

    lo = cells[icell].lo;
    hi = cells[icell].hi;

    ncreate = ncreate_values[icell];

    // if surfs in cell, use xcell for all created particle attempts

    if (cells[icell].nsurf)
      pflag = grid->point_outside_surfs(icell,xcell);

    for (int m = 0; m < ncreate; m++) {

      // generate random position X for new particle

      x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
      x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
      x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
      if (dimension == 2) x[2] = 0.0;

      // if surfs, check if random position is in flow region
      // if subcell, also check if in correct subcell
      // if not, attempt new insertion up to MAXATTEMPT times

      if (cells[icell].nsurf && pflag) {
        int nattempt = 0;
        while (nattempt < MAXATTEMPT) {
          if (grid->outside_surfs(icell,x,xcell)) {
            if (cells[icell].nsplit == 1) break;
            int splitcell = sinfo[cells[icell].isplit].icell;
            if (dimension == 2) subcell = update->split2d(splitcell,x);
            else subcell = update->split3d(splitcell,x);
            if (subcell == icell) break;
          }

          nattempt++;

          x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
          x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
          x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
          if (dimension == 2) x[2] = 0.0;
        }

        // particle insertion was unsuccessful

        if (nattempt >= MAXATTEMPT) continue;
      }

      // if region defined, skip if particle outside region

      if (region && !region->match(x)) continue;

      // insertion of particle at position X is accepted
      // calculate all other particle properties

      rn = random->uniform();

      isp = 0;
      while (cummulative[isp] < rn) isp++;
      ispecies = species[isp];

      if (speciesflag) {
        isp = species_variable(x) - 1;
        if (isp < 0 || isp >= nspecies) continue;
        ispecies = species[isp];
      }

      if (tempflag) {
        tempscale = temperature_variable(x);
        sqrttempscale = sqrt(tempscale);
      }

      vn = vscale[isp] * sqrttempscale * sqrt(-log(random->uniform()));
      vr = vscale[isp] * sqrttempscale * sqrt(-log(random->uniform()));
      theta1 = MY_2PI * random->uniform();
      theta2 = MY_2PI * random->uniform();

      if (velflag) {
        velocity_variable(x,vstream,vstream_variable);
        v[0] = vstream_variable[0] + vn*cos(theta1);
        v[1] = vstream_variable[1] + vr*cos(theta2);
        v[2] = vstream_variable[2] + vr*sin(theta2);
      } else {
        v[0] = vstream[0] + vn*cos(theta1);
        v[1] = vstream[1] + vr*cos(theta2);
        v[2] = vstream[2] + vr*sin(theta2);
      }

      erot = particle->erot(ispecies,temp_rot*tempscale,random);
      evib = particle->evib(ispecies,temp_vib*tempscale,random);

      id = MAXSMALLINT*random->uniform();

      particle->add_particle(id,ispecies,icell,x,v,erot,evib);

      if (nfix_update_custom)
        modify->update_custom(particle->nlocal-1,temp_thermal,
                             temp_rot,temp_vib,vstream);
    }
  }

  // clean up

  memory->destroy(ncreate_values);

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
   use particle position in svar variable to generate particle species
   first plug in particle x,y,z values into sxvar,syvar,szvar
------------------------------------------------------------------------- */

int CreateParticles::species_variable(double *x)
{
  if (sxstr) input->variable->internal_set(sxvar,x[0]);
  if (systr) input->variable->internal_set(syvar,x[1]);
  if (szstr) input->variable->internal_set(szvar,x[2]);

  double value = input->variable->compute_equal(svar);
  int isp = static_cast<int> (value);
  return isp;
}

/* ----------------------------------------------------------------------
   use grid cell center in dvar variable to generate density scale factor
   first plug in grid x,y,z values into dxvar,dyvar,dzvar
------------------------------------------------------------------------- */

double CreateParticles::density_variable(double *lo, double *hi)
{
  double center[3];
  center[0] = 0.5 * (lo[0]+hi[0]);
  center[1] = 0.5 * (lo[1]+hi[1]);
  center[2] = 0.5 * (lo[2]+hi[2]);

  if (dxstr) input->variable->internal_set(dxvar,center[0]);
  if (dystr) input->variable->internal_set(dyvar,center[1]);
  if (dzstr) input->variable->internal_set(dzvar,center[2]);

  double scale = input->variable->compute_equal(dvar);
  return scale;
}

/* ----------------------------------------------------------------------
   use grid cell center in tvar variable to generate temperature scale factor
   first plug in grid x,y,z values into txvar,tyvar,tzvar
------------------------------------------------------------------------- */

double CreateParticles::temperature_variable(double *x)
{

  if (txstr) input->variable->internal_set(txvar,x[0]);
  if (tystr) input->variable->internal_set(tyvar,x[1]);
  if (tzstr) input->variable->internal_set(tzvar,x[2]);

  double scale = input->variable->compute_equal(tvar);
  return scale;
}
/* ----------------------------------------------------------------------
   use particle position in vxvar,vyvar,vzvar variables to generate vel stream
   first plug in particle x,y,z values into vvarx,vvary,vvarz
------------------------------------------------------------------------- */

void CreateParticles::velocity_variable(double *x, double *vstream,
                                        double *vstream_variable)
{
  if (vstrx) input->variable->internal_set(vvarx,x[0]);
  if (vstry) input->variable->internal_set(vvary,x[1]);
  if (vstrz) input->variable->internal_set(vvarz,x[2]);

  if (vxstr) vstream_variable[0] = input->variable->compute_equal(vxvar);
  else vstream_variable[0] = vstream[0];
  if (vystr) vstream_variable[1] = input->variable->compute_equal(vyvar);
  else vstream_variable[1] = vstream[1];
  if (vzstr) vstream_variable[2] = input->variable->compute_equal(vzvar);
  else vstream_variable[2] = vstream[2];
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
  RanKnuth *random = new RandomPark(update->ranmaster->uniform());

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
