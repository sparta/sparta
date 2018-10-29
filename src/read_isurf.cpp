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
#include "string.h"
#include "stdlib.h"
#include "ctype.h"
#include "read_isurf.h"
#include "math_extra.h"
#include "surf.h"
#include "domain.h"
#include "grid.h"
#include "comm.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NEITHER,BAD,GOOD};
enum{NONE,CHECK,KEEP};
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files

#define MAXLINE 256
#define CHUNK 1024
#define EPSILON_NORM 1.0e-12
#define EPSILON_GRID 1.0e-3
#define BIG 1.0e20
#define DELTA 128           // must be 2 or greater 

// NOTE: allow reading 2nd set of isurfs into a different group region ??
// NOTE: check that all boundary point values are 0
// NOTE: option to write out surf once formed?

/* ---------------------------------------------------------------------- */

ReadISurf::ReadISurf(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);
  line = new char[MAXLINE];
  keyword = new char[MAXLINE];
  buffer = new char[CHUNK*MAXLINE];
}

/* ---------------------------------------------------------------------- */

ReadISurf::~ReadISurf()
{
  delete [] line;
  delete [] keyword;
  delete [] buffer;
}

/* ---------------------------------------------------------------------- */

void ReadISurf::command(int narg, char **arg)
{
  if (!grid->exist)
    error->all(FLERR,"Cannot read_isurf before grid is defined");
  if (!surf->implicit)
    error->all(FLERR,"Cannot read_isurf unless global surf implicit is set");
  if (particle->exist)
    error->all(FLERR,"Cannot read_isurf when particles exist");

  surf->exist = 1;
  dim = domain->dimension;

  if (narg != 2) error->all(FLERR,"Illegal read_isurf command");

  // read header info
  // verify that grid is uniform child cells within region

  if (me == 0) {
    if (screen) fprintf(screen,"Reading isurf file ...\n");
    open(arg[0]);
  }

  int igroup = grid->find_group(arg[1]);
  if (igroup < 0) error->all(FLERR,"Read_isurf group ID does not exist");

  int nx,ny,nz;
  //int flag = grid->check_uniform(igroup,nx,ny,nz);
  // NOTE: how to store offsets, some way to loop over grid cells in group

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  header();

  // read corner point grid values
  // assign 8 corner point values to each grid cell

  read_grid();

  // close file

  if (me == 0) {
    if (compressed) pclose(fp);
    else fclose(fp);
  }

  // convert grid point values to surfs
  // NOTE: this is where will call marching cubes?

  grid2surf();

  // extent of surfs
  // compute sizes of smallest surface elements

  /*
  double extent[3][2];
  extent[0][0] = extent[1][0] = extent[2][0] = BIG;
  extent[0][1] = extent[1][1] = extent[2][1] = -BIG;

  int m = npoint_old;
  for (int i = 0; i < npoint_new; i++) {
    extent[0][0] = MIN(extent[0][0],pts[m].x[0]);
    extent[0][1] = MAX(extent[0][1],pts[m].x[0]);
    extent[1][0] = MIN(extent[1][0],pts[m].x[1]);
    extent[1][1] = MAX(extent[1][1],pts[m].x[1]);
    extent[2][0] = MIN(extent[2][0],pts[m].x[2]);
    extent[2][1] = MAX(extent[2][1],pts[m].x[2]);
    m++;
  }

  double minlen,minarea;
  if (dim == 2) minlen = shortest_line();
  if (dim == 3) smallest_tri(minlen,minarea);

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  %g %g xlo xhi\n",extent[0][0],extent[0][1]);
      fprintf(screen,"  %g %g ylo yhi\n",extent[1][0],extent[1][1]);
      fprintf(screen,"  %g %g zlo zhi\n",extent[2][0],extent[2][1]);
      if (dim == 2)
	fprintf(screen,"  %g min line length\n",minlen);
      if (dim == 3) {
	fprintf(screen,"  %g min triangle edge length\n",minlen);
	fprintf(screen,"  %g min triangle area\n",minarea);
      }
    }
    if (logfile) {
      fprintf(logfile,"  %g %g xlo xhi\n",extent[0][0],extent[0][1]);
      fprintf(logfile,"  %g %g ylo yhi\n",extent[1][0],extent[1][1]);
      fprintf(logfile,"  %g %g zlo zhi\n",extent[2][0],extent[2][1]);
      if (dim == 2)
	fprintf(logfile,"  %g min line length\n",minlen);
      if (dim == 3) {
	fprintf(logfile,"  %g min triangle edge length\n",minlen);
	fprintf(logfile,"  %g min triangle area\n",minarea);
      }
    }
  }
  */

  // compute normals of new lines or triangles

  //if (dim == 2) surf->compute_line_normal(nline_old,nline_new);
  //else surf->compute_tri_normal(ntri_old,ntri_new);

  // error check on new points,lines,tris
  // all points must be inside or on surface of simulation box

  //surf->check_point_inside(npoint_old,npoint_new);

  // -----------------------
  // map surfs to grid cells
  // -----------------------

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // make list of surf elements I own
  // assign surfs to grid cells
  // error checks to flag bad surfs

  surf->setup_surf();

  grid->unset_neighbors();
  grid->remove_ghosts();
  grid->clear_surf();

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // error checks that can be done before surfs are mapped to grid cells

  if (dim == 2) {
    //surf->check_watertight_2d(npoint_old,nline_old);
  } else {
    //surf->check_watertight_3d(npoint_old,ntri_old);
  }

  MPI_Barrier(world);
  double time4 = MPI_Wtime();

  // map surfs to grid cells then error check
  // check done on per-grid-cell basis, too expensive to do globally

  grid->surf2grid(1);

  MPI_Barrier(world);
  double time5 = MPI_Wtime();

  // re-setup grid ghosts and neighbors

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->reset_neighbors();
  comm->reset_neighbors();

  MPI_Barrier(world);
  double time6 = MPI_Wtime();

  // flag cells and corners as OUTSIDE or INSIDE

  grid->set_inout();
  grid->type_check();

  MPI_Barrier(world);
  double time7 = MPI_Wtime();

  // stats

  double time_total = time6-time1;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  read/sort/check/surf2grid/ghost/inout/particle percent = "
              "%g %g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total,100.0*(time7-time6)/time_total);
    }
    if (logfile) {
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  read/sort/check/surf2grid/ghost/inout/particle percent = "
              "%g %g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total,100.0*(time7-time6)/time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   read free-format header of data file
   1st line and blank lines are skipped
   non-blank lines are checked for header keywords and leading value is read
   header ends with non-blank line containing no header keyword (or EOF)
   return line with non-blank line (or empty line if EOF)
------------------------------------------------------------------------- */

void ReadISurf::header()
{
  int n;
  char *ptr;

  // skip 1st line of file

  if (me == 0) {
    char *eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
  }

  //npoint_new = nline_new = ntri_new = 0;

  while (1) {

    // read a line and bcast length

    if (me == 0) {
      if (fgets(line,MAXLINE,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
    }
    MPI_Bcast(&n,1,MPI_INT,0,world);

    // if n = 0 then end-of-file so return with blank line

    if (n == 0) {
      line[0] = '\0';
      return;
    }

    // bcast line

    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // trim anything from '#' onward
    // if line is blank, continue

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

  }

  //if (npoint_new == 0) error->all(FLERR,"Surf file does not contain points");
}

/* ----------------------------------------------------------------------
   read/store all grid values
------------------------------------------------------------------------- */

void ReadISurf::read_grid()
{
  int i,m,nchunk;
  char *next,*buf;

  // read and broadcast one CHUNK of lines at a time

  /*

  int n = npoint_old;
  int nread = 0;
  
  while (nread < npoint_new) {
    if (npoint_new-nread > CHUNK) nchunk = CHUNK;
    else nchunk = npoint_new-nread;
    if (me == 0) {
      char *eof;
      m = 0;
      for (i = 0; i < nchunk; i++) {
	eof = fgets(&buffer[m],MAXLINE,fp);
	if (eof == NULL) error->one(FLERR,"Unexpected end of surf file");
	m += strlen(&buffer[m]);
      }
      m++;
    }
    MPI_Bcast(&m,1,MPI_INT,0,world);
    MPI_Bcast(buffer,m,MPI_CHAR,0,world);

    buf = buffer;
    next = strchr(buf,'\n');
    *next = '\0';
    int nwords = input->count_words(buf);
    *next = '\n';

    if (dim == 2 && nwords != 3)
      error->all(FLERR,"Incorrect point format in surf file");
    if (dim == 3 && nwords != 4)
      error->all(FLERR,"Incorrect point format in surf file");

    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');
      strtok(buf," \t\n\r\f");
      pts[n].x[0] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
      pts[n].x[1] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
      if (dim == 3) 
        pts[n].x[2] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
      else pts[n].x[2] = 0.0;
      n++;
      buf = next + 1;
    }

    nread += nchunk;
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  %d points\n",npoint_new);
    if (logfile) fprintf(logfile,"  %d points\n",npoint_new);
  }

  */
}

/* ----------------------------------------------------------------------
   create implicit surfs from grid point values
------------------------------------------------------------------------- */

void ReadISurf::grid2surf()
{
}

/* ----------------------------------------------------------------------
   return shortest line length
------------------------------------------------------------------------- */

double ReadISurf::shortest_line()
{
  /*
  double len = BIG;
  int m = nline_old;
  for (int i = 0; i < nline_new; i++) {
    len = MIN(len,surf->line_size(m));
    m++;
  }
  return len;
  */
  return 0.0;
}

/* ----------------------------------------------------------------------
   return shortest tri edge and smallest tri area
------------------------------------------------------------------------- */

void ReadISurf::smallest_tri(double &len, double &area)
{
  double lenone,areaone;

  /*
  len = area = BIG;
  int m = ntri_old;
  for (int i = 0; i < ntri_new; i++) {
    areaone = surf->tri_size(m,lenone);
    len = MIN(len,lenone);
    area = MIN(area,areaone);
    m++;
  }
  */
}

/* ----------------------------------------------------------------------
   proc 0 opens data file
   test if gzipped
------------------------------------------------------------------------- */

void ReadISurf::open(char *file)
{
  compressed = 0;
  char *suffix = file + strlen(file) - 3;
  if (suffix > file && strcmp(suffix,".gz") == 0) compressed = 1;
  if (!compressed) fp = fopen(file,"r");
  else {
#ifdef SPARTA_GZIP
    char gunzip[128];
    sprintf(gunzip,"gunzip -c %s",file);
    fp = popen(gunzip,"r");
#else
    error->one(FLERR,"Cannot open gzipped file");
#endif
  }

  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",file);
    error->one(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   grow surface data structures
------------------------------------------------------------------------- */

void ReadISurf::grow_surf()
{
  pts = (Surf::Point *) 
    memory->srealloc(pts,maxpoint*sizeof(Surf::Point),"surf:pts");
  lines = (Surf::Line *) 
    memory->srealloc(lines,maxline*sizeof(Surf::Line),"surf:lines");
  tris = (Surf::Tri *) 
    memory->srealloc(tris,maxtri*sizeof(Surf::Tri),"surf:tris");
}

