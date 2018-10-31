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
#include "write_surf.h"
#include "surf.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

#define MAXLINE 256

/* ---------------------------------------------------------------------- */

WriteSurf::WriteSurf(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

void WriteSurf::command(int narg, char **arg)
{
  if (!surf->exist)
    error->all(FLERR,"Cannot write surf when surfs do not exist");

  if (narg != 1) error->all(FLERR,"Illegal write_surf command");

  // write file

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  int me = comm->me;
  FILE *fp;

  if (me == 0) {
    if (screen) fprintf(screen,"Writing surf file ...\n");
    fp = fopen(arg[0],"w");
    if (!fp) {
      char str[128];
      sprintf(str,"Cannot open surface file %s",arg[0]);
      error->one(FLERR,str);
    }
  }

  if (me == 0) write_file(fp);

  // close file

  if (me == 0) fclose(fp);

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // stats

  double time_total = time2-time1;

  int nsurf;
  if (domain->dimension == 2) nsurf = surf->nline;
  else nsurf = surf->ntri;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  surf elements = %d\n",nsurf);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
    }

    if (logfile) {
      fprintf(logfile,"  surf elements = %d\n",nsurf);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   write surf file
   only called by proc 0
------------------------------------------------------------------------- */

void WriteSurf::write_file(FILE *fp)
{
  int dim = domain->dimension;

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  int nline = surf->nline;
  int ntri = surf->ntri;

  // header section

  fprintf(fp,"# Surface element file written by SPARTA\n\n");
  // NOTE POINT: how to do this
  //fprintf(fp,"%d points\n",npoint);
  if (dim == 2) fprintf(fp,"%d lines\n",nline);
  else fprintf(fp,"%d triangles\n",ntri);
  fprintf(fp,"\n");

  // points
  // NOTE POINT: how to do this?

  /*
  fprintf(fp,"Points\n\n");
  if (dim == 2) {
    for (int i = 0; i < npoint; i++)
      fprintf(fp,"%d %20.15g %20.15g\n",i+1,pts[i].x[0],pts[i].x[1]);
  } else {
    for (int i = 0; i < npoint; i++)
      fprintf(fp,"%d %20.15g %20.15g %20.15g\n",i+1,
	      pts[i].x[0],pts[i].x[1],pts[i].x[2]);
  }
  */

  // lines

  if (dim == 2) {
    fprintf(fp,"\nLines\n\n");
    for (int i = 0; i < nline; i++)
      fprintf(fp,"%d %d %d %d\n",i+1,lines[i].type,
	      lines[i].p1+1,lines[i].p2+1);
  }

  // triangles

  if (dim == 3) {
    fprintf(fp,"\nTriangles\n\n");
    for (int i = 0; i < ntri; i++)
      fprintf(fp,"%d %d %d %d %d\n",i+1,tris[i].type,
	      tris[i].p1+1,tris[i].p2+1,tris[i].p3+1);
  }
}
