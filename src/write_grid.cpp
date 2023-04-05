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
#include "write_grid.h"
#include "grid.h"
#include "domain.h"
#include "comm.h"
#include "hash3.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{IDS,GEOM};
#define MAXLINE 256

/* ---------------------------------------------------------------------- */

WriteGrid::WriteGrid(SPARTA *sparta) : Pointers(sparta)
{
  silent = 0;
}

/* ---------------------------------------------------------------------- */

void WriteGrid::command(int narg, char **arg)
{
  if (!grid->exist)
    error->all(FLERR,"Cannot write grid when grid is not defined");

  if (narg != 1) error->all(FLERR,"Illegal write_grid command");

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  // open file on proc 0

  int me = comm->me;
  if (me == 0) {
    if (screen && !silent) fprintf(screen,"Writing grid file ...\n");
    fp = fopen(arg[0],"w");
    if (!fp) {
      char str[128];
      sprintf(str,"Cannot open file %s",arg[0]);
      error->one(FLERR,str);
    }
  }

  // write file

  if (me == 0) header();
  write();

  // close file

  if (me == 0) fclose(fp);

  // stats

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  double time_total = time2-time1;

  if (comm->me == 0 && !silent) {
    if (screen) {
      fprintf(screen,"  grid cells = " BIGINT_FORMAT "\n",grid->ncell);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
    }

    if (logfile) {
      fprintf(logfile,"  grid cells = " BIGINT_FORMAT "\n",grid->ncell);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   write header of grid file
   only called by proc 0
------------------------------------------------------------------------- */

void WriteGrid::header()
{
  Grid::ParentLevel *plevels = grid->plevels;

  fprintf(fp,"# Grid file of cell IDs written by SPARTA\n\n");
  fprintf(fp,BIGINT_FORMAT " cells\n",grid->ncell);
  fprintf(fp,"%d levels\n",grid->maxlevel);
  for (int ilevel = 0; ilevel < grid->maxlevel; ilevel++)
    fprintf(fp,"%d %d %d level-%d\n",
            plevels[ilevel].nx,plevels[ilevel].ny,plevels[ilevel].nz,ilevel+1);
}

/* ----------------------------------------------------------------------
   write Cells section of grid file
   proc 0 pings other procs for info and writes entire file
------------------------------------------------------------------------- */

void WriteGrid::write()
{
  int i,tmp,nlines;
  MPI_Status status;
  MPI_Request request;

  int me = comm->me;
  int nprocs = comm->nprocs;

  // nme = # of cells this proc will contribute

  int nme = grid->nunsplitlocal + grid->nsplitlocal;

  // allocate memory for max of nme

  int nmax;
  MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);

  bigint *buf;
  memory->create(buf,nmax,"write_grid:buf");

  // pack ID each cell into buf, skipping sub cells

  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  int m = 0;
  for (i = 0; i < nglocal; i++) {
    if (cells[i].nsplit <= 0) continue;
    buf[m++] = cells[i].id;
  }

  // write out cell IDs from all procs

  if (me == 0) {
    fprintf(fp,"\nCells\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,nmax,MPI_SPARTA_BIGINT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&nlines);
      } else nlines = nme;

      for (i = 0; i < nlines; i++) {
        fprintf(fp,BIGINT_FORMAT "\n",buf[i]);
      }
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buf,nme,MPI_SPARTA_BIGINT,0,0,world);
  }

  // clean up

  memory->destroy(buf);
}
