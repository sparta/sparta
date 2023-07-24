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

  if (narg < 1) error->all(FLERR,"Illegal write_grid command");

  // optional args

  int iarg = 1;

  ncustom = 0;
  index_custom = NULL;
  type_custom = NULL;
  size_custom = NULL;
  
  while (iarg < narg) {
    if (strcmp(arg[iarg],"custom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid write_grid command");

      memory->grow(index_custom,ncustom+1,"writegrid:index_custom");
      memory->grow(type_custom,ncustom+1,"writegrid:type_custom");
      memory->grow(size_custom,ncustom+1,"writegridf:size_custom");

      int index = grid->find_custom(arg[iarg+1]);
      if (index < 0) error->all(FLERR,"Write_grid custom name does not exist");
      index_custom[ncustom] = index;
      type_custom[ncustom] = grid->etype[index];
      size_custom[ncustom] = grid->esize[index];
      ncustom++;

      iarg += 2;
    } else error->all(FLERR,"Invalid write_grid command");
  }
  
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

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  if (me == 0) header();
  write();

  // close file

  if (me == 0) fclose(fp);

  // clean up custom data

  if (ncustom) {
    memory->destroy(index_custom);
    memory->destroy(type_custom);
    memory->destroy(size_custom);
  }
  
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
  int i,ic,tmp,nlines;
  MPI_Status status,cstatus;
  MPI_Request request,crequest;

  int me = comm->me;
  int nprocs = comm->nprocs;

  // nme = # of cells this proc will contribute

  int nme = grid->nunsplitlocal + grid->nsplitlocal;

  // allocate memory for max of nme in idbuf and cbuf
  // nvalues_custom = # of custom values per grid cell

  int nmax;
  MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  
  bigint *idbuf;
  memory->create(idbuf,nmax,"write_grid:idbuf");

  nvalues_custom = 0;
  double **cbuf = NULL;

  if (ncustom) {
    for (ic = 0; ic < ncustom; ic++)
      if (size_custom[ic] == 0) nvalues_custom++;
      else nvalues_custom += size_custom[ic];
    memory->create(cbuf,nmax,nvalues_custom,"write_grid:cbuf");
  }
  
  // pack ID of each child cell into idbuf, skipping sub cells
  // pack custom values into cbuf
  
  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  int m = 0;
  for (i = 0; i < nglocal; i++) {
    if (cells[i].nsplit <= 0) continue;
    idbuf[m++] = cells[i].id;
  }

  if (ncustom) {
    int m = 0;
    for (i = 0; i < nglocal; i++) {
      if (cells[i].nsplit <= 0) continue;
      pack_custom(i,cbuf[m]);
      m++;
    }
  }
  
  // write out cell IDs and custom values from all procs

  if (me == 0) {
    fprintf(fp,"\nCells\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(idbuf,nmax,MPI_SPARTA_BIGINT,iproc,0,world,&request);
	if (ncustom) MPI_Irecv(&cbuf[0][0],nmax*nvalues_custom,MPI_DOUBLE,
			       iproc,0,world,&crequest);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
	if (ncustom) MPI_Wait(&crequest,&cstatus);
        MPI_Get_count(&status,MPI_SPARTA_BIGINT,&nlines);
      } else nlines = nme;

      for (i = 0; i < nlines; i++) {
        fprintf(fp,BIGINT_FORMAT,idbuf[i]);
	if (ncustom) write_custom(cbuf[i]);
	fprintf(fp,"\n");
      }
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(idbuf,nme,MPI_SPARTA_BIGINT,0,0,world);
    if (ncustom) {
      if (nme) MPI_Rsend(&cbuf[0][0],nme,MPI_DOUBLE,0,0,world);
      else MPI_Rsend(NULL,nme,MPI_DOUBLE,0,0,world);
    }
  }

  // clean up

  memory->destroy(idbuf);
  if (ncustom) memory->destroy(cbuf);
}

/* ----------------------------------------------------------------------
   pack user-specified custom values for Ith grid cell into vec
---------------------------------------------------------------------- */

void WriteGrid::pack_custom(int i, double *vec)
{
  int m = 0;
  
  for (int ic = 0; ic < ncustom; ic++) {
    if (type_custom[ic] == 0) {
      if (size_custom[ic] == 0) {
	int *ivector = grid->eivec[grid->ewhich[index_custom[ic]]];
	vec[m++] = ivector[i];
      } else {
	int **iarray = grid->eiarray[grid->ewhich[index_custom[ic]]];
	for (int j = 0; j < size_custom[ic]; j++)
	  vec[m++] = iarray[i][j];
      }
    } else {
      if (size_custom[ic] == 0) {
	double *dvector = grid->edvec[grid->ewhich[index_custom[ic]]];
	vec[m++] = dvector[i];
      } else {
	double **darray = grid->edarray[grid->ewhich[index_custom[ic]]];
	for (int j = 0; j < size_custom[ic]; j++)
	  vec[m++] = darray[i][j];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   write user-specified custom values from vec to output file
   vec has ncustom_values for a single grid cell
---------------------------------------------------------------------- */

void WriteGrid::write_custom(double *vec)
{
  int m = 0;
  
  for (int ic = 0; ic < ncustom; ic++) {
    if (type_custom[ic] == 0) {
      if (size_custom[ic] == 0) {
	fprintf(fp," %d",(int) vec[m++]);
      } else {
	for (int j = 0; j < size_custom[ic]; j++)
	  fprintf(fp," %d",(int) vec[m++]);
      }
    } else {
      if (size_custom[ic] == 0) {
	fprintf(fp," %g",vec[m++]);
      } else {
	for (int j = 0; j < size_custom[ic]; j++)
	  fprintf(fp," %g",vec[m++]);
      }
    }
  }
}
