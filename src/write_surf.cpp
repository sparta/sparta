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

  if (narg < 1) error->all(FLERR,"Illegal write_surf command");

  // optional args

  int pflag = 1;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[1],"points") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal write_surf command");
      if (strcmp(arg[2],"yes") == 0) pflag = 1;
      else if (strcmp(arg[2],"no") == 0) pflag = 0;
      else error->all(FLERR,"Illegal write_surf command");
      iarg += 2;
    } else error->all(FLERR,"Illegal write_surf command");
  }

  // open file

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  FILE *fp;

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Writing surf file ...\n");
    fp = fopen(arg[0],"w");
    if (!fp) {
      char str[128];
      sprintf(str,"Cannot open surface file %s",arg[0]);
      error->one(FLERR,str);
    }
  }

  // write file

  write_file(fp,pflag);

  // close file

  if (comm->me == 0) fclose(fp);

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

/* ---------------------------------------------------------------------- */

void WriteSurf::write_file(FILE *fp, int pflag)
{
  me = comm->me;
  nprocs = comm->nprocs;

  if (!surf->distributed && pflag) {
    if (me == 0) write_file_all_points(fp);
  } else if (!surf->distributed && !pflag) {
    if (me == 0) write_file_all_nopoints(fp);
  } else if (pflag) {
    write_file_distributed_points(fp);
  } else if (!pflag) {
    write_file_distributed_nopoints(fp);
  }
}

/* ----------------------------------------------------------------------
   write surf file for explicit non-distributed elements
   only called by proc 0
   write Point section with duplicate points
------------------------------------------------------------------------- */

void WriteSurf::write_file_all_points(FILE *fp)
{
  int dim = domain->dimension;

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  int nsurf = surf->nlocal;

  // header section

  fprintf(fp,"# Surface element file written by SPARTA\n\n");
  if (dim == 2) {
    bigint n = (bigint) 2 * nsurf;
    fprintf(fp,BIGINT_FORMAT " points\n",n);
    fprintf(fp,"%d lines\n",nsurf);
  } else if (dim == 3) {
    bigint n = (bigint) 3 * nsurf;
    fprintf(fp,BIGINT_FORMAT " points\n",n);
    fprintf(fp,"%d triangles\n",nsurf);
  }

  fprintf(fp,"\nPoints\n\n");

  if (dim == 2) {
    for (int i = 0; i < nsurf; i++) {
      fprintf(fp,BIGINT_FORMAT " %20.15g %20.15g\n",
              (bigint) 2*i+1,lines[i].p1[0],lines[i].p1[1]);
      fprintf(fp,BIGINT_FORMAT " %20.15g %20.15g\n",
              (bigint) 2*i+2,lines[i].p2[0],lines[i].p2[1]);
    }
  } else {
    for (int i = 0; i < nsurf; i++) {
      fprintf(fp,BIGINT_FORMAT " %20.15g %20.15g %20.15g\n",
              (bigint) 3*i+1,tris[i].p1[0],tris[i].p1[1],tris[i].p1[2]);
      fprintf(fp,BIGINT_FORMAT " %20.15g %20.15g %20.15g\n",
              (bigint) 3*i+2,tris[i].p2[0],tris[i].p2[1],tris[i].p2[2]);
      fprintf(fp,BIGINT_FORMAT " %20.15g %20.15g %20.15g\n",
              (bigint) 3*i+3,tris[i].p3[0],tris[i].p3[1],tris[i].p3[2]);
    }
  }

  // lines

  if (dim == 2) {
    fprintf(fp,"\nLines\n\n");
    for (int i = 0; i < nsurf; i++)
      fprintf(fp,SURFINT_FORMAT " %d " BIGINT_FORMAT " " BIGINT_FORMAT "\n",
              lines[i].id,lines[i].type,
              (bigint) 2*i+1, (bigint) 2*i+2);
  }
  
  // triangles

  if (dim == 3) {
    fprintf(fp,"\nTriangles\n\n");
    for (int i = 0; i < nsurf; i++)
      fprintf(fp,SURFINT_FORMAT " %d " BIGINT_FORMAT " " 
              BIGINT_FORMAT " " BIGINT_FORMAT "\n",
              tris[i].id,tris[i].type,
              (bigint) 3*i+1, (bigint) 3*i+2, (bigint) 3*i+3);
  }
}

/* ----------------------------------------------------------------------
   write surf file for explicit non-distributed elements
   only called by proc 0
   no Point section, include point coords with Lines/Triangles
------------------------------------------------------------------------- */

void WriteSurf::write_file_all_nopoints(FILE *fp)
{
  int dim = domain->dimension;

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  int nsurf = surf->nlocal;

  // header section

  fprintf(fp,"# Surface element file written by SPARTA\n\n");
  if (dim == 2) fprintf(fp,"%d lines\n",nsurf);
  else fprintf(fp,"%d triangles\n",nsurf);

  // lines

  if (dim == 2) {
    fprintf(fp,"\nLines\n\n");
    for (int i = 0; i < nsurf; i++)
      fprintf(fp,SURFINT_FORMAT " %d %20.15g %20.15g %20.15g %20.15g\n",
              lines[i].id,lines[i].type,
              lines[i].p1[0],lines[i].p1[1],
              lines[i].p2[0],lines[i].p2[1]);
  }
  
  // triangles

  if (dim == 3) {
    fprintf(fp,"\nTriangles\n\n");
    for (int i = 0; i < nsurf; i++)
      fprintf(fp,SURFINT_FORMAT " %d %20.15g %20.15g %20.15g "
              "%20.15g %20.15g %20.15g %20.15g %20.15g %20.15g\n",
              tris[i].id,tris[i].type,
              tris[i].p1[0],tris[i].p1[1],tris[i].p1[2],
              tris[i].p2[0],tris[i].p2[1],tris[i].p2[2],
              tris[i].p3[0],tris[i].p3[1],tris[i].p3[2]);
  }
}

/* ----------------------------------------------------------------------
   write surf file for explicit or implicit distributed elements
   all procs participate
   write Point section with duplicate points
------------------------------------------------------------------------- */

void WriteSurf::write_file_distributed_points(FILE *fp)
{
  int dim = domain->dimension;
  bigint nsurf = surf->nsurf;

  // nmine = my implicit or explicit surf count

  int nmine;
  if (surf->implicit) nmine = surf->nlocal;
  else nmine = surf->nown;

  // check for overflow

  bigint ndouble,nchar;

  if (dim == 2) {
    ndouble = (bigint) 2*nmine * 2;
    nchar = (bigint) nmine * sizeof(SurfIDType);
  } else {
    ndouble = (bigint) 3*nmine * 3;
    nchar = (bigint) nmine * sizeof(SurfIDType);
  }

  printf("WRITESURF %d\n",nmine);


  if (ndouble > MAXSMALLINT || nchar > MAXSMALLINT) 
    error->one(FLERR,"Too much distributed data to communicate");

  // header section

  if (me == 0) {
    fprintf(fp,"# Surface element file written by SPARTA\n\n");
    if (dim == 2) {
      fprintf(fp,BIGINT_FORMAT " points\n",(bigint) 2*nsurf);
      fprintf(fp,BIGINT_FORMAT " lines\n",nsurf);
    } else {
      fprintf(fp,BIGINT_FORMAT " points\n",(bigint) 3*nsurf);
      fprintf(fp,BIGINT_FORMAT " triangles\n",nsurf);
    }
  }

  // each proc contributes explicit or implicit distributed points
  // nmine = element count, npoint = point count
  // nper = size of ID + type
  // lines_mine/tris_mine = ptr in Surf to elements

  int npoint,nper;
  Surf::Line *lines_mine;
  Surf::Tri *tris_mine;

  if (surf->implicit) {
    if (dim == 2) {
      npoint = 2*nmine;
      lines_mine = surf->lines;
    } else {
      npoint = 3*nmine;
      tris_mine = surf->tris;
    }
  } else {
    if (dim == 2) {
      npoint = 2*nmine;
      lines_mine = surf->mylines;
    } else {
      npoint = 3*nmine;
      tris_mine = surf->mytris;
    }
  }
  nper = sizeof(SurfIDType);

  // max_size_point = largest point buffer needed by any proc

  int max_size_point;
  MPI_Allreduce(&npoint,&max_size_point,1,MPI_INT,MPI_MAX,world);

  // pbuf = local buf for my explicit or implicit points

  double *pbuf;
  memory->create(pbuf,max_size_point*dim,"writesurf:pbuf");

  // pack my points into buf

  int m = 0;

  if (dim == 2) {
    for (int i = 0; i < nmine; i++) {
      pbuf[m++] = lines_mine[i].p1[0];
      pbuf[m++] = lines_mine[i].p1[1];
      pbuf[m++] = lines_mine[i].p2[0];
      pbuf[m++] = lines_mine[i].p2[1];
    }
  } else {
    for (int i = 0; i < nmine; i++) {
      pbuf[m++] = tris_mine[i].p1[0];
      pbuf[m++] = tris_mine[i].p1[1];
      pbuf[m++] = tris_mine[i].p1[2];
      pbuf[m++] = tris_mine[i].p2[0];
      pbuf[m++] = tris_mine[i].p2[1];
      pbuf[m++] = tris_mine[i].p2[2];
      pbuf[m++] = tris_mine[i].p3[0];
      pbuf[m++] = tris_mine[i].p3[1];
      pbuf[m++] = tris_mine[i].p3[2];
    }
  }

  // if proc 0: ping each proc, receive its point data, write data to file
  // else: wait for ping from proc 0, send my data to proc 0

  int tmp,recv_size,ncount;
  MPI_Status status;
  MPI_Request request;

  if (me == 0) fprintf(fp,"\nPoints\n\n");

  if (me == 0) {
    bigint index = 0;
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(pbuf,max_size_point*dim,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&recv_size);
      } else recv_size = npoint*dim;
      
      ncount = recv_size/dim;
      m = 0;

      if (dim == 2) {
	for (int i = 0; i < ncount; i++) {
	  index++;
	  fprintf(fp,BIGINT_FORMAT " %20.15g %20.15g\n",
                  index,pbuf[m],pbuf[m+1]);
          m += 2;
	}
      } else {
	for (int i = 0; i < ncount; i++) {
	  index++;
	  fprintf(fp,BIGINT_FORMAT " %20.15g %20.15g %20.15g\n",
                  index,pbuf[m],pbuf[m+1],pbuf[m+2]);
          m += 3;
	}
      }
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(pbuf,npoint*dim,MPI_DOUBLE,0,0,world);
  }
  
  memory->sfree(pbuf);

  // max_size_surf = largest surf buffer needed by any proc

  int max_size_surf;
  MPI_Allreduce(&nmine,&max_size_surf,1,MPI_INT,MPI_MAX,world);

  // sbuf = local buf for my surface IDs and types

  SurfIDType *sbuf = (SurfIDType *) 
    memory->smalloc(max_size_surf*nper,"writesurf:sbuf");

  // pack my line/tri IDs/types into sbuf

  m = 0;
  if (dim == 2) {
    for (int i = 0; i < nmine; i++) {
      sbuf[i].id = lines_mine[i].id;
      sbuf[i].type = lines_mine[i].type;
    }
  } else {
    for (int i = 0; i < nmine; i++) {
      sbuf[i].id = tris_mine[i].id;
      sbuf[i].type = tris_mine[i].type;
    }
  }

  // if proc 0: ping each proc, receive its point data, write data to file
  // else: wait for ping from proc 0, send my data to proc 0

  if (me == 0) {
    if (dim == 2) fprintf(fp,"\nLines\n\n");
    else fprintf(fp,"\nTriangles\n\n");
  }

  if (me == 0) {
    bigint index = 0;
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(sbuf,max_size_surf*nper,MPI_CHAR,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_CHAR,&recv_size);
      } else recv_size = nmine*nper;
      
      ncount = recv_size/nper;
      if (dim == 2) {
	for (int i = 0; i < ncount; i++) {
	  fprintf(fp,SURFINT_FORMAT " %d " BIGINT_FORMAT " " BIGINT_FORMAT "\n",
		  sbuf[i].id,sbuf[i].type,index+1,index+2);
	  index += 2;
	}
      } else {
	for (int i = 0; i < ncount; i++) {
	  fprintf(fp,SURFINT_FORMAT " %d " BIGINT_FORMAT " " BIGINT_FORMAT " " 
		  BIGINT_FORMAT "\n",
		  sbuf[i].id,sbuf[i].type,index+1,index+2,index+3);
	  index += 3;
	}
      }
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(sbuf,nmine*nper,MPI_CHAR,0,0,world);
  }
  
  memory->sfree(sbuf);
}

/* ----------------------------------------------------------------------
   write surf file for explicit or implicit distributed elements
   all procs participate
   no Point section, include point coords with Lines/Triangles
------------------------------------------------------------------------- */

void WriteSurf::write_file_distributed_nopoints(FILE *fp)
{
  int dim = domain->dimension;
  bigint nsurf = surf->nsurf;

  // nmine = my implicit or explicit surf count

  int nmine;
  if (surf->implicit) nmine = surf->nlocal;
  else nmine = surf->nown;

  // check for overflow

  bigint nchar;

  if (dim == 2) nchar = (bigint) nmine * sizeof(Surf::Line);
  else nchar = (bigint) nmine * sizeof(Surf::Tri);

  if (nchar > MAXSMALLINT) 
    error->one(FLERR,"Too much distributed data to communicate");

  // header section

  if (me == 0) {
    fprintf(fp,"# Surface element file written by SPARTA\n\n");
    if (dim == 2)
      fprintf(fp,BIGINT_FORMAT " lines\n",nsurf);
    else if (dim == 3)
      fprintf(fp,BIGINT_FORMAT " triangles\n",nsurf);
  }

  // each proc contributes explicit or implicit distributed surfs
  // nmine = element count
  // lines_mine/tris_mine = ptr in Surf to elements

  int nper;
  Surf::Line *lines_mine;
  Surf::Tri *tris_mine;

  if (dim == 2) nper = sizeof(Surf::Line);
  else nper = sizeof(Surf::Tri);

  if (surf->implicit) {
    if (dim == 2) lines_mine = surf->lines;
    else tris_mine = surf->tris;
  } else {
    if (dim == 2) lines_mine = surf->mylines;
    else tris_mine = surf->mytris;
  }

  // max_size = largest buffer needed by any proc

  int max_size;
  MPI_Allreduce(&nmine,&max_size,1,MPI_INT,MPI_MAX,world);

  // lines/tris = local buf for my explicit or implicit elements
  // pack my surfs into buf

  Surf::Line *lines;
  Surf::Tri *tris;
  char *buf;

  if (dim == 2) {
    lines = (Surf::Line *) memory->smalloc(max_size*nper,"writesurf:lines");
    memcpy(lines,lines_mine,nmine*nper);
    buf = (char *) lines;
  } else  {
    tris = (Surf::Tri *) memory->smalloc(max_size*nper,"writesurf:tris");
    memcpy(tris,tris_mine,nmine*nper);
    buf = (char *) tris;
  }

  // if proc 0: ping each proc, receive its surf data, write data to file
  // else: wait for ping from proc 0, send my data to proc 0

  int tmp,recv_size,ncount;
  MPI_Status status;
  MPI_Request request;

  if (me == 0) {
    if (dim == 2) fprintf(fp,"\nLines\n\n");
    else fprintf(fp,"\nTriangles\n\n");
  }

  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,max_size*nper,MPI_CHAR,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_CHAR,&recv_size);
      } else recv_size = nmine*nper;
      
      ncount = recv_size/nper;
      if (dim == 2) {
	lines = (Surf::Line *) buf;
	for (int i = 0; i < ncount; i++)
	  fprintf(fp,SURFINT_FORMAT " %d %20.15g %20.15g %20.15g %20.15g\n",
		  lines[i].id,lines[i].type,
		  lines[i].p1[0],lines[i].p1[1],
		  lines[i].p2[0],lines[i].p2[1]);
      } else {
	tris = (Surf::Tri *) buf;
	for (int i = 0; i < ncount; i++)
	  fprintf(fp,SURFINT_FORMAT " %d %20.15g %20.15g %20.15g "
                  "%20.15g %20.15g %20.15g %20.15g %20.15g %20.15g\n",
		  tris[i].id,tris[i].type,
		  tris[i].p1[0],tris[i].p1[1],tris[i].p1[2],
		  tris[i].p2[0],tris[i].p2[1],tris[i].p2[2],
		  tris[i].p3[0],tris[i].p3[1],tris[i].p3[2]);
      }
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buf,nmine*nper,MPI_CHAR,0,0,world);
  }
  
  memory->sfree(buf);
}
