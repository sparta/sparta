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

#include "mpi.h"
#include "spatype.h"
#include "string.h"
#include "write_isurf.h"
#include "update.h"
#include "domain.h"
#include "grid.h"
#include "surf.h"
#include "comm.h"
#include "modify.h"
#include "fix.h"
#include "fix_ablate.h"
#include "input.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{INT,DOUBLE};

/* ---------------------------------------------------------------------- */

WriteISurf::WriteISurf(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

void WriteISurf::command(int narg, char **arg)
{
  MPI_Comm_rank(world,&me);
  dim = domain->dimension;

  if (!surf->exist || surf->implicit == 0)
    error->all(FLERR,"Cannot write_isurf when implicit surfs do not exist");

  if (narg < 6) error->all(FLERR,"Illegal write_isurf command");

  ggroup = grid->find_group(arg[0]);
  if (ggroup < 0) error->all(FLERR,"Write_isurf grid group ID does not exist");
  groupbit = grid->bitmask[ggroup];

  nx = input->inumeric(FLERR,arg[1]);
  ny = input->inumeric(FLERR,arg[2]);
  nz = input->inumeric(FLERR,arg[3]);

  if (dim == 2 && nz != 1) error->all(FLERR,"Invalid write_isurf command");

  // if filename contains a "*", replace with current timestep

  char *ptr;
  int n = strlen(arg[4]) + 16;
  char *file = new char[n];

  if ((ptr = strchr(arg[4],'*'))) {
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s",arg[4],update->ntimestep,ptr+1);
  } else strcpy(file,arg[4]);

  // ablation fix ID

  char *ablateID = arg[5];
  int ifix = modify->find_fix(ablateID);
  if (ifix < 0)
    error->all(FLERR,"Fix ID for write_isurf does not exist");
  if (strcmp(modify->fix[ifix]->style,"ablate") != 0)
    error->all(FLERR,"Fix for write_surf is not a fix ablate");
  ablate = (FixAblate *) modify->fix[ifix];

  // check that group and Nx,Ny,Nz match FixAblate

  if (ggroup != ablate->igroup)
    error->all(FLERR,"Write_isurf group does not match fix ablate group");

  if (nx != ablate->nx || ny != ablate->ny || nz != ablate->nz)
    error->all(FLERR,"Write_isurf Nxyz does not match fix ablate Nxyz");

  // process optional command line args

  precision = DOUBLE;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"precision") == 0)  {
      if (iarg+2 > narg) error->all(FLERR,"Invalid write_isurf command");
      if (strcmp(arg[iarg+1],"int") == 0) precision = INT;
      else if (strcmp(arg[iarg+1],"double") == 0) precision = DOUBLE;
      else error->all(FLERR,"Invalid write_isurf command");
      iarg += 2;
    } else error->all(FLERR,"Invalid write_isurf command");
  }

  // collect all corner point data into one big vector

  if (me == 0 && screen) fprintf(screen,"Writing isurf file ...\n");

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  dbuf = dbufall = NULL;

  collect_values();

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // write file

  FILE *fp;

  if (me == 0) {
    fp = fopen(file,"wb");
    if (!fp) {
      char str[128];
      sprintf(str,"Cannot open grid corner point file %s",file);
      error->one(FLERR,str);
    }
  }

  if (me == 0) write_file(fp);

  // close file

  if (me == 0) fclose(fp);

  memory->destroy(dbuf);
  memory->destroy(dbufall);

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // stats

  double time_total = time3-time1;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  corner points = %d\n",ncorner);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  collect/write percent = %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total);
    }

    if (logfile) {
      fprintf(logfile,"  corner points = %d\n",ncorner);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  collect/write percent = %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   collect corner point values into one big array
   each proc owns copy, use MPI_Allreduce to sum
------------------------------------------------------------------------- */

void WriteISurf::collect_values()
{
  double *cornerlo = ablate->cornerlo;
  double *xyzsize = ablate->xyzsize;

  bigint bncorner = bigint (nx+1) * (ny+1);
  if (dim == 3) bncorner *= (nz+1);
  if (bncorner > MAXSMALLINT) error->all(FLERR,"Write_isurf grid is too large");
  ncorner = bncorner;

  memory->create(dbuf,ncorner,"write_isurf:dbuf");
  memory->create(dbufall,ncorner,"write_isurf:dbufall");
  for (int i = 0; i < ncorner; i++) dbuf[i] = 0.0;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  int ix,iy,iz,index;
  double **array_grid = ablate->array_grid;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    // ix,iy,iz = cell index (1 to Nxyz) within array of grid cells

    ix =
      static_cast<int> ((cells[icell].lo[0]-cornerlo[0]) / xyzsize[0] + 0.5) + 1;
    iy =
      static_cast<int> ((cells[icell].lo[1]-cornerlo[1]) / xyzsize[1] + 0.5) + 1;
    iz =
      static_cast<int> ((cells[icell].lo[2]-cornerlo[2]) / xyzsize[2] + 0.5) + 1;

    // index = corner point index, 0 to (Nx+1)*(Ny+1)*(Nz+1) - 1
    // x varies fastest, then y, z slowest
    // this is for lower/left/bottom corner point of icell's 4/8 points

    index = (iz-1)*(ny+1)*(nx+1) + (iy-1)*(nx+1) + ix-1;

    // always copy lower/left/bottom corner point
    // copy other corner points if on upper boundary of Nx,Ny,Nz

    dbuf[index] = array_grid[icell][0];
    if (ix == nx) dbuf[index+1] = array_grid[icell][1];
    if (iy == ny) dbuf[index+nx+1] = array_grid[icell][2];
    if (ix == nx && iy == ny) dbuf[index+nx+2] = array_grid[icell][3];

    if (iz == nz && dim == 3) {
      index += (ny+1)*(nx+1);
      dbuf[index] = array_grid[icell][4];
      if (ix == nx) dbuf[index+1] = array_grid[icell][5];
      if (iy == ny) dbuf[index+nx+1] = array_grid[icell][6];
      if (ix == nx && iy == ny) dbuf[index+nx+2] = array_grid[icell][7];
    }
  }

  // MPI_Allreduce to sum dbuf across all procs
  // so that proc 0 has a copy to write to file

  MPI_Allreduce(dbuf,dbufall,ncorner,MPI_DOUBLE,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   write grid corner point file
   only called by proc 0
------------------------------------------------------------------------- */

void WriteISurf::write_file(FILE *fp)
{
  int nxyz[3];
  nxyz[0] = nx + 1;
  nxyz[1] = ny + 1;
  nxyz[2] = nz + 1;
  fwrite(nxyz,sizeof(int),dim,fp);

  if (precision == INT) {
    uint8_t *ibuf8;
    memory->create(ibuf8,ncorner,"write_isurf:ibuf8");
    for (int i = 0; i < ncorner; i++)
      ibuf8[i] = static_cast<int> (dbufall[i]);
    fwrite(ibuf8,sizeof(uint8_t),ncorner,fp);
    memory->destroy(ibuf8);
  }

  if (precision == DOUBLE) fwrite(dbufall,sizeof(double),ncorner,fp);
}
