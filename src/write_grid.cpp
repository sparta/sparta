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

#include "mpi.h"
#include "spatype.h"
#include "string.h"
#include "write_grid.h"
#include "grid.h"
#include "domain.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#ifdef SPARTA_MAP
#include <map>
#else
#include <tr1/unordered_map>
#endif

using namespace SPARTA_NS;

enum{PARENT,GEOM};
#define MAXLINE 256

/* ---------------------------------------------------------------------- */

WriteGrid::WriteGrid(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

void WriteGrid::command(int narg, char **arg)
{
  if (!grid->exist)
    error->all(FLERR,"Cannot write grid when grid is not defined");

  if (narg != 2) error->all(FLERR,"Illegal write_grid command");

  int mode;
  if (strcmp(arg[0],"parent") == 0) mode = PARENT;
  else if (strcmp(arg[0],"geom") == 0) mode = GEOM;
  else error->all(FLERR,"Illegal write_grid command");

  // write file, create parent cells and then child cells

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  int me = comm->me;
  if (me == 0) {
    if (screen) fprintf(screen,"Writing grid file ...\n");
    fp = fopen(arg[1],"w");
    if (!fp) {
      char str[128];
      sprintf(str,"Cannot open file %s",arg[1]);
      error->one(FLERR,str);
    }
  }

  // write Parents section

  if (mode == PARENT) {
    if (me == 0) header_parents();
    if (me == 0) write_parents();
  } else if (mode == GEOM) {
    if (me == 0) header_geometry();
    write_geometry();
  }

  // close file

  if (me == 0) fclose(fp);

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // stats

  double time_total = time2-time1;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  parent cells = %d\n",grid->nparent);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
    }

    if (logfile) {
      fprintf(logfile,"  parent cells = %d\n",grid->nparent);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   write header of parent grid file
   only called by proc 0
------------------------------------------------------------------------- */

void WriteGrid::header_parents()
{
  fprintf(fp,"# Parent grid file written by SPARTA\n\n");
  fprintf(fp,"%d nparents\n",grid->nparent);
}

/* ----------------------------------------------------------------------
   write Parents section of grid file
   only called by proc 0
------------------------------------------------------------------------- */

void WriteGrid::write_parents()
{
  char str[32];

  // fill hash with parent IDs if necessary

  if (!grid->hashfilled) {

#ifdef SPARTA_MAP
    std::map<cellint,int> *hash = grid->hash;
#else
    std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

    hash->clear();

    Grid::ParentCell *pcells = grid->pcells;
    int nparent = grid->nparent;

    for (int icell = 0; icell < nparent; icell++)
      (*hash)[pcells[icell].id] = -(icell+1);
  }

  fprintf(fp,"\nParents\n\n");

  Grid::ParentCell *pcells = grid->pcells;
  int nparent = grid->nparent;

  // one parent cell per line

  for (int i = 0; i < nparent; i++) {
    grid->id_num2str(pcells[i].id,str);
    fprintf(fp,"%d %s %d %d %d\n",i+1,str,
            pcells[i].nx,pcells[i].ny,pcells[i].nz);
  }

  // clear hash if filled it

  if (!grid->hashfilled) grid->hash->clear();
}

/* ----------------------------------------------------------------------
   write header of geometry grid file
   only called by proc 0
------------------------------------------------------------------------- */

void WriteGrid::header_geometry()
{
  fprintf(fp,"# Geometry grid file written by SPARTA\n\n");

  if (domain->dimension == 2) 
    fprintf(fp,BIGINT_FORMAT " points\n",4*grid->ncell);
  else fprintf(fp,BIGINT_FORMAT " points\n",8*grid->ncell);
  fprintf(fp,BIGINT_FORMAT " cells\n",grid->ncell);
}

/* ----------------------------------------------------------------------
   write Points,Cells sections of geometry grid file
   proc 0 pings other procs for info and writes entire file
------------------------------------------------------------------------- */

void WriteGrid::write_geometry()
{
  int i,tmp,nlines;
  double *buf;
  MPI_Status status;
  MPI_Request request;

  int me = comm->me;
  int nprocs = comm->nprocs;
  int dimension = domain->dimension;

  // nme = # of points this proc will contribute
  // NOTE: need to worry about 8*ncells overflowing

  int nme = grid->nunsplitlocal + grid->nsplitlocal;
  if (dimension == 2) nme *= 4;
  else nme *= 8;

  // nmax = max # of points on any proc

  int nmax;
  MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);

  // allocate memory for max # of points

  double **pt;
  memory->create(pt,nmax,3,"write_grid:pt");
  if (pt) buf = &pt[0][0];
  else buf = NULL;

  // pack corner points of each cell into pt, skipping sub cells

  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  int m = 0;
  for (i = 0; i < nglocal; i++) {
    if (cells[i].nsplit <= 0) continue;

    pt[m][0] = cells[i].lo[0];
    pt[m][1] = cells[i].lo[1];
    pt[m][2] = cells[i].lo[2];
    m++;

    pt[m][0] = cells[i].hi[0];
    pt[m][1] = cells[i].lo[1];
    pt[m][2] = cells[i].lo[2];
    m++;

    pt[m][0] = cells[i].hi[0];
    pt[m][1] = cells[i].hi[1];
    pt[m][2] = cells[i].lo[2];
    m++;

    pt[m][0] = cells[i].lo[0];
    pt[m][1] = cells[i].hi[1];
    pt[m][2] = cells[i].lo[2];
    m++;

    if (dimension == 2) continue;

    pt[m][0] = cells[i].lo[0];
    pt[m][1] = cells[i].lo[1];
    pt[m][2] = cells[i].hi[2];
    m++;

    pt[m][0] = cells[i].hi[0];
    pt[m][1] = cells[i].lo[1];
    pt[m][2] = cells[i].hi[2];
    m++;

    pt[m][0] = cells[i].hi[0];
    pt[m][1] = cells[i].hi[1];
    pt[m][2] = cells[i].hi[2];
    m++;

    pt[m][0] = cells[i].lo[0];
    pt[m][1] = cells[i].hi[1];
    pt[m][2] = cells[i].hi[2];
    m++;
  }

  // write out points from all procs

  printf("Ncells %d %d\n",me,nme);

  if (me == 0) {
    fprintf(fp,"\nPoints\n\n");
    bigint index = 0;
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,nmax*3,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&nlines);
        nlines /= 3;
      } else nlines = nme;
      
      for (i = 0; i < nlines; i++) {
        index++;
        fprintf(fp,BIGINT_FORMAT " %g %g %g\n",
                index,pt[i][0],pt[i][1],pt[i][2]);
      }
    }
    
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buf,nme*3,MPI_DOUBLE,0,0,world);
  }

  // dellocate point memory

  memory->destroy(pt);

  // proc 0 write Cells section
  // trivial monotonically increasing indexing b/c Points are not unique

  if (me != 0) return;

  bigint ncell = grid->ncell;
  bigint np = 0;

  fprintf(fp,"\nCells\n\n");

  if (dimension == 2) {
    for (bigint n = 0; n < ncell; n++) {
      fprintf(fp,BIGINT_FORMAT " " BIGINT_FORMAT " " BIGINT_FORMAT 
              " " BIGINT_FORMAT " " BIGINT_FORMAT "\n",
              n+1,np+1,np+2,np+3,np+4);
      np += 4;
    }
  } else {
    for (bigint n = 0; n < ncell; n++) {
      fprintf(fp,BIGINT_FORMAT " " BIGINT_FORMAT " " BIGINT_FORMAT 
              " " BIGINT_FORMAT " " BIGINT_FORMAT
              " " BIGINT_FORMAT " " BIGINT_FORMAT
              " " BIGINT_FORMAT " " BIGINT_FORMAT "\n",
              n+1,np+1,np+2,np+3,np+4,np+5,np+6,np+7,np+8);
      np += 8;
    }
  }
}
