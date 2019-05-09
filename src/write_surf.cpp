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

  // 3 kinds of surfaces

  if (!surf->distributed) {
    if (me == 0) write_file_all(fp);
  } else if (surf->distributed && !surf->implcit) {
    write_file_distributed(fp);
  } else if (surf->implicit) {
    write_file_implicit(fp);
  }

  // close file

  if (me == 0) fclose(fp);

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // stats

  double time_total = time2-time1;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  surf elements = " BIGINT_FORMAT "\n",surf->nsurf);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
    }

    if (logfile) {
      fprintf(logfile,"  surf elements = " BIGINT_FORMAT "\n",surf->nsurf);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   write surf file for explicit non-distributed elements
   only called by proc 0
   write points with duplicates
   each line writes its 2 end points
   each triangle writes its 3 corner points
------------------------------------------------------------------------- */

void WriteSurf::write_file_all(FILE *fp)
{
  int dim = domain->dimension;

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  // for explicit & non-distributed, Nsurf will fit in 32-bit integer

  int nsurf = surf->nlocal;

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

/* ----------------------------------------------------------------------
   write surf file for explicit distributed elements
   all procs participate
   write points with duplicates
   each line writes its 2 end points
   each triangle writes its 3 corner points
------------------------------------------------------------------------- */

void WriteSurf::write_file_distributed(FILE *fp)
{
  int dim = domain->dimension;

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  // for explicit & non-distributed, Nsurf will fit in 32-bit integer

  bigint nsurf = surf->nsurf;

  // header section

  fprintf(fp,"# Surface element file written by SPARTA\n\n");
  if (dim == 2) {
    fprintf(fp,BIGINT_FORMAT " points\n",2*nsurf);
    fprintf(fp,BIGINT_FORMAT " lines\n",nsurf);
  } else if (dim == 3) {
    fprintf(fp,BIGINT_FORMAT " points\n",3*nsurf);
    fprintf(fp,BIGINT_FORMAT " triangles\n",nsurf);
  }
  fprintf(fp,"\n");

  // communication buffer for point info
  // max_size = largest buffer needed by any proc
  // NOTE: check for overflow

  int nown = surf->nown;
  int send_size;
  if (dim = 2) send_size = 4*nown;
  else send_size = 9*nown;
  int max_size;
  MPI_Allreduce(&send_size,&max_size,1,MPI_INT,MPI_MAX,world);

  double *buf;
  memory->create(buf,max_size,"write_surf:buf");

  // pack my points into buf

  int m = 0;

  if (dim == 2) {
    for (int i = 0; i < nown; i++) {
      buf[m++] = surf->mylines[i].p1[0];
      buf[m++] = surf->mylines[i].p1[1];
      buf[m++] = surf->mylines[i].p2[0];
      buf[m++] = surf->mylines[i].p2[1];
    }
  } else {
    for (int i = 0; i < nown; i++) {
      buf[m++] = surf->mylines[i].p1[0];
      buf[m++] = surf->mylines[i].p1[1];
      buf[m++] = surf->mylines[i].p1[2];
      buf[m++] = surf->mylines[i].p2[0];
      buf[m++] = surf->mylines[i].p2[1];
      buf[m++] = surf->mylines[i].p2[2];
      buf[m++] = surf->mylines[i].p3[0];
      buf[m++] = surf->mylines[i].p3[1];
      buf[m++] = surf->mylines[i].p3[2];
    }
  }

  // if proc 0: ping each proc, receive its data, write data to file
  // else: wait for ping from proc 0, send my data to proc 0

  int tmp,recv_size,n;
  MPI_Status status;
  MPI_Request request;

  fprintf(fp,"Points\n\n");

  if (me == 0) {
    bigint index = 0;

    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,max_size,MPI_CHAR,me+iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_CHAR,&recv_size);
      } else recv_size = send_size;
      
      if (dim == 2) {
        n = recv_size / 4;
        m = 0;
        for (int i = 0; i < n; i++) {
          fprintf(fp,"%d %20.15g %20.15g\n",index+1,buf[m+0],buf[m+1]);
          fprintf(fp,"%d %20.15g %20.15g\n",index+2,buf[m+2],buf[m+3]);
          m += 4;
          index += 2;
        }
      } else {
        n = recv_size / 9;
        m = 0;
        for (int i = 0; i < nsurf; i++) {
          fprintf(fp,"%d %20.15g %20.15g %20.15g\n",
                  index+1,buf[m+0],buf[m+1],buf[m+2]);
          fprintf(fp,"%d %20.15g %20.15g %20.15g\n",
                  index+2,buf[m+3],buf[m+4],buf[m+5]);
          fprintf(fp,"%d %20.15g %20.15g %20.15g\n",
                  index+3,buf[m+6],buf[m+7],buf[m+8]);
          m += 9;
          index += 2;
        }
      }
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buf,send_size,MPI_CHAR,0,0,world);
  }

  memory->destroy(buf);

  // communication buffer for surface element type info
  // max_size = largest buffer needed by any proc
  // NOTE: check for overflow
  // NOTE: need surf IDs as well to keep them the same
  //       are IDs used when read them in?

  int nown = surf->nown;
  MPI_Allreduce(&nown,&max_size,1,MPI_INT,MPI_MAX,world);

  double *ibuf;
  memory->create(ibuf,max_size,"write_surf:ibuf");

  // pack my surf types into buf

  if (dim == 2) {
    for (int i = 0; i < nown; i++)
      ibuf[m++] = surf->mylines[i].type;
  } else {
    for (int i = 0; i < nown; i++)
      ibuf[m++] = surf->mytris[i].type;
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
