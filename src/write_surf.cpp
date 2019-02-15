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

  int nsurf = surf->nsurf;

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
   write points with duplicates
   each line writes its 2 end points
   each triangle writes its 3 corner points
------------------------------------------------------------------------- */

void WriteSurf::write_file(FILE *fp)
{
  int dim = domain->dimension;

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  int nsurf = surf->nsurf;

  // header section

  fprintf(fp,"# Surface element file written by SPARTA\n\n");
  if (dim == 2) {
    fprintf(fp,"%d points\n",2*nsurf);
    fprintf(fp,"%d lines\n",nsurf);
  } else if (dim == 3) {
    fprintf(fp,"%d points\n",3*nsurf);
    fprintf(fp,"%d triangles\n",nsurf);
  }
  fprintf(fp,"\n");

  fprintf(fp,"Points\n\n");

  if (dim == 2) {
    for (int i = 0; i < nsurf; i++) {
      fprintf(fp,"%d %20.15g %20.15g\n",2*i+1,lines[i].p1[0],lines[i].p1[1]);
      fprintf(fp,"%d %20.15g %20.15g\n",2*i+2,lines[i].p2[0],lines[i].p2[1]);
    }
  } else {
    for (int i = 0; i < nsurf; i++) {
      fprintf(fp,"%d %20.15g %20.15g %20.15g\n",
              3*i+1,tris[i].p1[0],tris[i].p1[1],tris[i].p1[2]);
      fprintf(fp,"%d %20.15g %20.15g %20.15g\n",
              3*i+2,tris[i].p2[0],tris[i].p2[1],tris[i].p2[2]);
      fprintf(fp,"%d %20.15g %20.15g %20.15g\n",
              3*i+3,tris[i].p3[0],tris[i].p3[1],tris[i].p3[2]);
    }
  }

  // lines

  if (dim == 2) {
    fprintf(fp,"\nLines\n\n");
    for (int i = 0; i < nsurf; i++)
      fprintf(fp,"%d %d %d %d\n",i+1,lines[i].type,2*i+1,2*i+2);
  }

  // triangles

  if (dim == 3) {
    fprintf(fp,"\nTriangles\n\n");
    for (int i = 0; i < nsurf; i++)
      fprintf(fp,"%d %d %d %d %d\n",i+1,tris[i].type,3*i+1,3*i+2,3*i+3);
  }
}
