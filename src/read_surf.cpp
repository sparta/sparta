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
#include "read_surf.h"
#include "math_extra.h"
#include "surf.h"
#include "domain.h"
#include "grid.h"
#include "comm.h"
#include "geometry.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

#ifdef SPARTA_MAP
#include <map>
#elif SPARTA_UNORDERED_MAP
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

using namespace SPARTA_NS;
using namespace MathConst;

enum{NEITHER,BAD,GOOD};

#define MAXLINE 256
#define CHUNK 1024
#define EPSILON_NORM 1.0e-12
#define EPSILON_GRID 1.0e-3
#define BIG 1.0e20
#define DELTA 128           // must be 2 or greater 

/* ---------------------------------------------------------------------- */

ReadSurf::ReadSurf(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);
  line = new char[MAXLINE];
  keyword = new char[MAXLINE];
  buffer = new char[CHUNK*MAXLINE];
}

/* ---------------------------------------------------------------------- */

ReadSurf::~ReadSurf()
{
  delete [] line;
  delete [] keyword;
  delete [] buffer;
}

/* ---------------------------------------------------------------------- */

void ReadSurf::command(int narg, char **arg)
{
  if (!grid->exist) 
    error->all(FLERR,"Cannot read_surf before grid is defined");
  if (!grid->exist_ghost)
    error->all(FLERR,"Cannot read_surf before grid ghost cells are defined");
  if (particle->exist) 
    error->all(FLERR,"Cannot read_surf after particles are defined");

  surf->exist = 1;
  dimension = domain->dimension;

  if (narg < 1) error->all(FLERR,"Illegal read_surf command");

  // read header info

  if (me == 0) {
    if (screen) fprintf(screen,"Reading surf file ...\n");
    open(arg[0]);
  }

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  header();

  // extend pts,lines,tris data structures

  pts = surf->pts;
  lines = surf->lines;
  tris = surf->tris;

  npoint_old = surf->npoint;
  nline_old = surf->nline;
  ntri_old = surf->ntri;

  maxpoint = npoint_old + npoint_new;
  maxline = nline_old + nline_new;
  maxtri = ntri_old + ntri_new;

  pts = (Surf::Point *) 
    memory->srealloc(pts,maxpoint*sizeof(Surf::Point),"surf:pts");
  lines = (Surf::Line *) 
    memory->srealloc(lines,maxline*sizeof(Surf::Line),"surf:lines");
  tris = (Surf::Tri *) 
    memory->srealloc(tris,maxtri*sizeof(Surf::Tri),"surf:tris");

  // read and store Points and Lines/Tris sections

  parse_keyword(1);
  if (strcmp(keyword,"Points") != 0)
    error->all(FLERR,
	       "Read_surf did not find points section of surf file");
  read_points();

  parse_keyword(0);
  if (dimension == 2) {
    if (strcmp(keyword,"Lines") != 0)
      error->all(FLERR,
	       "Read_surf did not find lines section of surf file");
    read_lines();
  } else {
    if (strcmp(keyword,"Triangles") != 0)
    error->all(FLERR,
	       "Read_surf did not find triangles section of surf file");
    read_tris();
  }

  // close file

  if (me == 0) {
    if (compressed) pclose(fp);
    else fclose(fp);
  }

  // apply optional keywords for geometric transformations
  // store optional keywords for group and type information

  origin[0] = origin[1] = origin[2] = 0.0;
  int grouparg = 0;
  int typeadd = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"origin") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Invalid read_surf command");
      double ox = atof(arg[iarg+1]);
      double oy = atof(arg[iarg+2]);
      double oz = atof(arg[iarg+3]);
      if (dimension == 2 && oz != 0.0) 
	error->all(FLERR,"Invalid read_surf geometry transformation "
		   "for 2d simulation");
      origin[0] = ox;
      origin[1] = oy;
      origin[2] = oz;
      iarg += 4;
    } else if (strcmp(arg[iarg],"trans") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Invalid read_surf command");
      double dx = atof(arg[iarg+1]);
      double dy = atof(arg[iarg+2]);
      double dz = atof(arg[iarg+3]);
      if (dimension == 2 && dz != 0.0) 
	error->all(FLERR,"Invalid read_surf geometry transformation "
		   "for 2d simulation");
      origin[0] += dx;
      origin[1] += dy;
      origin[2] += dz;
      translate(dx,dy,dz);
      iarg += 4;
    } else if (strcmp(arg[iarg],"atrans") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Invalid read_surf command");
      double ax = atof(arg[iarg+1]);
      double ay = atof(arg[iarg+2]);
      double az = atof(arg[iarg+3]);
      if (dimension == 2 && az != 0.0) 
	error->all(FLERR,"Invalid read_surf geometry transformation "
		   "for 2d simulation");
      double dx = ax - origin[0];
      double dy = ay - origin[1];
      double dz = az - origin[2];
      origin[0] = ax;
      origin[1] = ay;
      origin[2] = az;
      translate(dx,dy,dz);
      iarg += 4;
    } else if (strcmp(arg[iarg],"ftrans") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Invalid read_surf command");
      double fx = atof(arg[iarg+1]);
      double fy = atof(arg[iarg+2]);
      double fz = atof(arg[iarg+3]);
      if (dimension == 2 && fz != 0.5) 
	error->all(FLERR,"Invalid read_surf geometry transformation "
		   "for 2d simulation");
      double ax = domain->boxlo[0] + fx*domain->xprd;
      double ay = domain->boxlo[1] + fy*domain->yprd;
      double az;
      if (dimension == 3) az = domain->boxlo[2] + fz*domain->zprd;
      else az = 0.0;
      double dx = ax - origin[0];
      double dy = ay - origin[1];
      double dz = az - origin[2];
      origin[0] = ax;
      origin[1] = ay;
      origin[2] = az;
      translate(dx,dy,dz);
      iarg += 4;
    } else if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Invalid read_surf command");
      double sx = atof(arg[iarg+1]);
      double sy = atof(arg[iarg+2]);
      double sz = atof(arg[iarg+3]);
      if (dimension == 2 && sz != 1.0) 
	error->all(FLERR,"Invalid read_surf geometry transformation "
		   "for 2d simulation");
      scale(sx,sy,sz);
      iarg += 4;
    } else if (strcmp(arg[iarg],"rotate") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Invalid read_surf command");
      double theta = atof(arg[iarg+1]);
      double rx = atof(arg[iarg+2]);
      double ry = atof(arg[iarg+3]);
      double rz = atof(arg[iarg+4]);
      if (dimension == 2 && (rx != 0.0 || ry != 0.0 || rz != 1.0))
	error->all(FLERR,"Invalid read_surf geometry transformation "
		   "for 2d simulation");
      if (rx == 0.0 && ry == 0.0 && rz == 0.0)
	error->all(FLERR,"Invalid read_surf geometry transformation "
		   "for 2d simulation");
      rotate(theta,rx,ry,rz);
      iarg += 5;
    } else if (strcmp(arg[iarg],"invert") == 0) {
      invert();
      iarg += 1;
    } else if (strcmp(arg[iarg],"clip") == 0) {
      double frac = 0.0;
      if (iarg+1 < narg) {
        char c = arg[iarg+1][0];
        if (isdigit(c) || c == '-' || c == '+' || c == '.') {
          frac = atof(arg[iarg+1]);
          if (frac < 0.0 || frac >= 0.5) 
            error->all(FLERR,"Invalid read_surf command");
          iarg++;
        }
      }
      if (frac > 0.0) push_points_to_boundary(frac);
      if (dimension == 2) clip2d();
      else clip3d();
      iarg++;

    } else if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid read_surf command");
      grouparg = iarg+1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"typeadd") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid read_surf command");
      typeadd = input->inumeric(FLERR,arg[iarg+1]);
      if (typeadd < 0) error->all(FLERR,"Invalid read_surf command");
      iarg += 2;

    } else error->all(FLERR,"Invalid read_surf command");
  }

  // if specified, apply group and typeadd keywords
  // these reset per-element mask/type info

  if (grouparg) {
    int igroup = surf->find_group(arg[grouparg]);
    if (igroup < 0) igroup = surf->add_group(arg[grouparg]);
    int groupbit = surf->bitmask[igroup];
    int n = npoint_old + npoint_new;
    if (dimension == 2)
      for (int i = npoint_old; i < n; i++) lines[i].mask |= groupbit;
    else
      for (int i = npoint_old; i < n; i++) tris[i].mask |= groupbit;
  }

  if (typeadd) {
    int n = npoint_old + npoint_new;
    if (dimension == 2)
      for (int i = npoint_old; i < n; i++) lines[i].type += typeadd;
    else 
      for (int i = npoint_old; i < n; i++) tris[i].type += typeadd;
  }

  // update Surf data structures

  surf->pts = pts;
  surf->lines = lines;
  surf->tris = tris;

  surf->npoint = npoint_old + npoint_new;
  surf->nline = nline_old + nline_new;
  surf->ntri = ntri_old + ntri_new;

  // extent of surf after geometric transformations
  // compute sizes of smallest surface elements

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
  if (dimension == 2) minlen = shortest_line();
  if (dimension == 3) smallest_tri(minlen,minarea);

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  %g %g xlo xhi\n",extent[0][0],extent[0][1]);
      fprintf(screen,"  %g %g ylo yhi\n",extent[1][0],extent[1][1]);
      fprintf(screen,"  %g %g zlo zhi\n",extent[2][0],extent[2][1]);
      if (dimension == 2)
	fprintf(screen,"  %g min line length\n",minlen);
      if (dimension == 3) {
	fprintf(screen,"  %g min triangle edge length\n",minlen);
	fprintf(screen,"  %g min triangle area\n",minarea);
      }
    }
    if (logfile) {
      fprintf(logfile,"  %g %g xlo xhi\n",extent[0][0],extent[0][1]);
      fprintf(logfile,"  %g %g ylo yhi\n",extent[1][0],extent[1][1]);
      fprintf(logfile,"  %g %g zlo zhi\n",extent[2][0],extent[2][1]);
      if (dimension == 2)
	fprintf(logfile,"  %g min line length\n",minlen);
      if (dimension == 3) {
	fprintf(logfile,"  %g min triangle edge length\n",minlen);
	fprintf(logfile,"  %g min triangle area\n",minarea);
      }
    }
  }

  // compute normals of new lines or triangles

  if (dimension == 2) surf->compute_line_normal(nline_old,nline_new);
  if (dimension == 3) surf->compute_tri_normal(ntri_old,ntri_new);

  // error check on new points,lines,tris
  // all points must be inside or on surface of simulation box

  check_point_inside();

  // make list of surf elements I own
  // assign surfs to grid cells
  // more error checks to flag bad surfs
  // re-setup ghost and grid

  surf->setup_surf();

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  grid->unset_neighbors();
  grid->remove_ghosts();
  grid->clear_surf();

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  grid->surf2grid(1);

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  if (dimension == 2) {
    check_watertight_2d();
    check_neighbor_norm_2d();
    check_point_near_surf_2d();
  } else {
    check_watertight_3d();
    check_neighbor_norm_3d();
    check_point_near_surf_3d();
  }

  MPI_Barrier(world);
  double time4 = MPI_Wtime();

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->reset_neighbors();
  comm->reset_neighbors();

  MPI_Barrier(world);
  double time5 = MPI_Wtime();

  grid->set_inout();
  grid->type_check();

  // DEBUG
  //grid->debug();

  MPI_Barrier(world);
  double time6 = MPI_Wtime();

  double time_total = time6-time1;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  read/surf2grid/error/ghost/inout percent = "
              "%g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total);
    }
    if (logfile) {
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  read/surf2grid/error/ghost/inout percent = "
              "%g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total);
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

void ReadSurf::header()
{
  int n;
  char *ptr;

  // skip 1st line of file

  if (me == 0) {
    char *eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
  }

  npoint_new = nline_new = ntri_new = 0;

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

    if (ptr = strchr(line,'#')) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    // search line for header keyword and set corresponding variable

    if (strstr(line,"points")) sscanf(line,"%d",&npoint_new);
    else if (strstr(line,"lines")) {
      if (dimension == 3) 
	error->all(FLERR,"Surf file cannot contain lines for 3d simulation");
      sscanf(line,"%d",&nline_new);
    } else if (strstr(line,"triangles")) {
      if (dimension == 2) 
	error->all(FLERR,
		   "Surf file cannot contain triangles for 2d simulation");
      sscanf(line,"%d",&ntri_new);
    } else break;
  }

  if (npoint_new == 0) error->all(FLERR,"Surf file does not contain points");
  if (dimension == 2 && nline_new == 0) 
    error->all(FLERR,"Surf file does not contain lines");
  if (dimension == 3 && ntri_new == 0) 
    error->all(FLERR,"Surf file does not contain triangles");
}

/* ----------------------------------------------------------------------
   read/store all points
------------------------------------------------------------------------- */

void ReadSurf::read_points()
{
  int i,m,nchunk;
  char *next,*buf;

  // read and broadcast one CHUNK of lines at a time

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
    int nwords = count_words(buf);
    *next = '\n';

    if (dimension == 2 && nwords != 3)
      error->all(FLERR,"Incorrect point format in surf file");
    if (dimension == 3 && nwords != 4)
      error->all(FLERR,"Incorrect point format in surf file");

    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');
      strtok(buf," \t\n\r\f");
      pts[n].x[0] = atof(strtok(NULL," \t\n\r\f"));
      pts[n].x[1] = atof(strtok(NULL," \t\n\r\f"));
      if (dimension == 3) pts[n].x[2] = atof(strtok(NULL," \t\n\r\f"));
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
}

/* ----------------------------------------------------------------------
   read/store all lines
   alter p1,p2 indices to point to newest set of stored points
------------------------------------------------------------------------- */

void ReadSurf::read_lines()
{
  int i,m,nchunk,type,p1,p2;
  char *next,*buf;

  // read and broadcast one CHUNK of lines at a time

  int n = nline_old;
  int nread = 0;

  while (nread < nline_new) {
    if (nline_new-nread > CHUNK) nchunk = CHUNK;
    else nchunk = nline_new-nread;
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
    int nwords = count_words(buf);
    *next = '\n';

    // allow for optional type in each line element

    if (nwords != 3 && nwords != 4)
      error->all(FLERR,"Incorrect line format in surf file");
    int typeflag = 0;
    if (nwords == 4) typeflag = 1;

    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');
      strtok(buf," \t\n\r\f");
      if (typeflag) type = atoi(strtok(NULL," \t\n\r\f"));
      else type = 1;
      p1 = atoi(strtok(NULL," \t\n\r\f"));
      p2 = atoi(strtok(NULL," \t\n\r\f"));
      if (p1 < 1 || p1 > npoint_new || p2 < 1 || p2 > npoint_new || p1 == p2)
	error->all(FLERR,"Invalid point index in line");
      lines[n].type = type;
      lines[n].mask = 1;
      lines[n].isc = lines[n].isr = -1;
      lines[n].p1 = p1-1 + npoint_old;
      lines[n].p2 = p2-1 + npoint_old;
      n++;
      buf = next + 1;
    }

    nread += nchunk;
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  %d lines\n",nline_new);
    if (logfile) fprintf(logfile,"  %d lines\n",nline_new);
  }
}

/* ----------------------------------------------------------------------
   read/store all triangles
   alter p1,p2,p3 indices to point to newest set of stored points
------------------------------------------------------------------------- */

void ReadSurf::read_tris()
{
  int i,m,nchunk,type,p1,p2,p3;
  char *next,*buf;

  // read and broadcast one CHUNK of triangles at a time

  int n = ntri_old;
  int nread = 0;

  while (nread < ntri_new) {
    if (ntri_new-nread > CHUNK) nchunk = CHUNK;
    else nchunk = ntri_new-nread;
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
    int nwords = count_words(buf);
    *next = '\n';

    // allow for optional type in each tri element

    if (nwords != 4 && nwords != 5)
      error->all(FLERR,"Incorrect line format in surf file");
    int typeflag = 0;
    if (nwords == 5) typeflag = 1;

    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');
      strtok(buf," \t\n\r\f");
      if (typeflag) type = atoi(strtok(NULL," \t\n\r\f"));
      else type = 1;
      p1 = atoi(strtok(NULL," \t\n\r\f"));
      p2 = atoi(strtok(NULL," \t\n\r\f"));
      p3 = atoi(strtok(NULL," \t\n\r\f"));
      if (p1 < 1 || p1 > npoint_new || p2 < 1 || p2 > npoint_new || 
	  p3 < 1 || p3 > npoint_new || p1 == p2 || p2 == p3)
	error->all(FLERR,"Invalid point index in triangle");
      tris[n].type = type;
      tris[n].mask = 1;
      tris[n].isc = tris[n].isr = -1;
      tris[n].p1 = p1-1 + npoint_old;
      tris[n].p2 = p2-1 + npoint_old;
      tris[n].p3 = p3-1 + npoint_old;

      n++;
      buf = next + 1;
    }

    nread += nchunk;
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  %d triangles\n",ntri_new);
    if (logfile) fprintf(logfile,"  %d triangles\n",ntri_new);
  }
}

/* ----------------------------------------------------------------------
   translate new vertices by (dx,dy,dz)
   for 2d, dz will be 0.0
------------------------------------------------------------------------- */

void ReadSurf::translate(double dx, double dy, double dz)
{
  int m = npoint_old;
  for (int i = 0; i < npoint_new; i++) {
    pts[m].x[0] += dx;
    pts[m].x[1] += dy;
    pts[m].x[2] += dz;
    m++;
  }
}

/* ----------------------------------------------------------------------
   scale new vertices by (sx,sy,sz) around origin
   for 2d, do not reset x[2] to avoid epsilon change
------------------------------------------------------------------------- */

void ReadSurf::scale(double sx, double sy, double sz)
{
  int m = npoint_old;
  for (int i = 0; i < npoint_new; i++) {
    pts[m].x[0] = sx*(pts[m].x[0]-origin[0]) + origin[0];
    pts[m].x[1] = sy*(pts[m].x[1]-origin[1]) + origin[1];
    if (dimension == 3) pts[m].x[2] = sz*(pts[m].x[2]-origin[2]) + origin[2];
    m++;
  }
}

/* ----------------------------------------------------------------------
   rotate new vertices around origin
   for 2d, do not reset x[2] to avoid epsilon change
------------------------------------------------------------------------- */

void ReadSurf::rotate(double theta, double rx, double ry, double rz)
{
  double r[3],q[4],d[3],dnew[3];
  double rotmat[3][3];

  theta *= MY_PI/180.0;

  r[0] = rx; r[1] = ry; r[2] = rz;
  MathExtra::norm3(r);
  MathExtra::axisangle_to_quat(r,theta,q);
  MathExtra::quat_to_mat(q,rotmat);

  int m = npoint_old;
  for (int i = 0; i < npoint_new; i++) {
    d[0] = pts[m].x[0] - origin[0];
    d[1] = pts[m].x[1] - origin[1];
    d[2] = pts[m].x[2] - origin[2];
    MathExtra::matvec(rotmat,d,dnew);
    pts[m].x[0] = dnew[0] + origin[0];
    pts[m].x[1] = dnew[1] + origin[1];
    if (dimension == 3) pts[m].x[2] = dnew[2] + origin[2];
    m++;
  }
}

/* ----------------------------------------------------------------------
   invert new vertex ordering within each line or tri
   this flips direction of surface normal
------------------------------------------------------------------------- */

void ReadSurf::invert()
{
  int tmp;

  if (dimension == 2) {
    int m = nline_old;
    for (int i = 0; i < nline_new; i++) {
      tmp = lines[m].p1;
      lines[m].p1 = lines[m].p2;
      lines[m].p2 = tmp;
      m++;
    }
  }

  if (dimension == 3) {
    int m = ntri_old;
    for (int i = 0; i < ntri_new; i++) {
      tmp = tris[m].p2;
      tris[m].p2 = tris[m].p3;
      tris[m].p3 = tmp;
      m++;
    }
  }
}

/* ----------------------------------------------------------------------
   clip all lines so fit inside simulation box
   may discard some lines and their points completely
   new points and lines may be added which touch box surface
   condense data structures by removing deleted points & lines
------------------------------------------------------------------------- */

void ReadSurf::clip2d()
{
  int i,m,n,dim,side,flag;
  double value,param;
  double *x,*x1,*x2;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  for (int iface = 0; iface < 4; iface++) {
    dim = iface / 2;
    side = iface % 2;
    if (side == 0) value = boxlo[dim];
    else value = boxhi[dim];

    // check if line straddles clipping edge
    // straddle = pts on different sides (not on) clipping edge
    // create new point on edge and split line into two lines
    // reallocate pts and lines if needed and reset x1,x2

    m = nline_old;
    n = nline_new;

    for (i = 0; i < n; i++) {
      x1 = pts[lines[m].p1].x;
      x2 = pts[lines[m].p2].x;
      flag = 0;
      if (x1[dim] < value && x2[dim] > value) flag = 1;
      if (x1[dim] > value && x2[dim] < value) flag = 1;
      if (flag) {
	if (npoint_new == maxpoint) {
	  maxpoint += DELTA;
	  pts = (Surf::Point *) 
	    memory->srealloc(pts,maxpoint*sizeof(Surf::Point),"surf:pts");
	  x1 = pts[lines[m].p1].x;
	  x2 = pts[lines[m].p2].x;
	}
	if (nline_new == maxline) {
	  maxline += DELTA;
	  lines = (Surf::Line *) 
	    memory->srealloc(lines,maxline*sizeof(Surf::Line),"surf:lines");
	}

	param = (value-x1[dim]) / (x2[dim]-x1[dim]);
	pts[npoint_new].x[0] = x1[0] + param*(x2[0]-x1[0]);
	pts[npoint_new].x[1] = x1[1] + param*(x2[1]-x1[1]);
	pts[npoint_new].x[2] = 0.0;
	pts[npoint_new].x[dim] = value;
	npoint_new++;

	memcpy(&lines[nline_new],&lines[m],sizeof(Surf::Line));
	lines[m].p2 = npoint_new-1;
	lines[nline_new].p1 = npoint_new-1;
	nline_new++;
      }
      m++;
    }

    // project all points outside clipping edge to the edge

    m = npoint_old;
    for (i = 0; i < npoint_new; i++) {
      x = pts[m].x;
      if (side == 0 && x[dim] < value) x[dim] = value;
      if (side == 1 && x[dim] > value) x[dim] = value;
      m++;
    }

    // ptflag[I] = # of lines that include point I
    
    int *ptflag;
    memory->create(ptflag,npoint_new+npoint_old,"readsurf:ptflag");
    m = npoint_old;
    for (i = 0; i < npoint_new; i++) ptflag[m++] = 0;

    m = nline_old;
    for (i = 0; i < nline_new; i++) {
      ptflag[lines[m].p1]++;
      ptflag[lines[m].p2]++;
      m++;
    }

    // lineflag[I] = 1 if line I can be removed b/c lies on clip edge
    // also decrement ptflag for pts in removed line
    
    int *lineflag;
    memory->create(lineflag,nline_new+nline_old,"readsurf:lineflag");
    m = nline_old;
    for (i = 0; i < nline_new; i++) lineflag[m++] = 0;

    m = nline_old;
    for (i = 0; i < nline_new; i++) {
      if (pts[lines[m].p1].x[dim] == value && 
	  pts[lines[m].p2].x[dim] == value) {
	lineflag[m] = 1;
	ptflag[lines[m].p1]--;
	ptflag[lines[m].p2]--;
      }
      m++;
    }

    // condense pts/lines to remove deleted points and lines
    // when delete points, set ptflag[I] = new index of point I for kept points
    // renumber point indices in lines using altered ptflag
    // reset npoint_new,nline_new

    n = m = npoint_old;
    for (i = 0; i < npoint_new; i++) {
      if (ptflag[m]) {
	memcpy(&pts[n],&pts[m],sizeof(Surf::Point));
	ptflag[m] = n;
	n++;
      }
      m++;
    }
    npoint_new = n - npoint_old;
    
    n = m = nline_old;
    for (i = 0; i < nline_new; i++) {
      if (!lineflag[m]) {
	memcpy(&lines[n],&lines[m],sizeof(Surf::Line));
	lines[n].p1 = ptflag[lines[n].p1];
	lines[n].p2 = ptflag[lines[n].p2];
	n++;
      }
      m++;
    }
    nline_new = n - nline_old;

    // clean up

    memory->destroy(ptflag);
    memory->destroy(lineflag);
  }

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  clipped to %d points\n",npoint_new);
      fprintf(screen,"  clipped to %d lines\n",nline_new);
    }
    if (logfile) {
      fprintf(logfile,"  clipped to %d points\n",npoint_new);
      fprintf(logfile,"  clipped to %d lines\n",nline_new);
    }
  }
}

/* ----------------------------------------------------------------------
   clip all tris so fit inside simulation box
   may discard some tris and their points completely
   new points and tris may be added which touch box surface
   condense data structures by removing deleted points & tris
------------------------------------------------------------------------- */

void ReadSurf::clip3d()
{
  int i,m,n,dim,side,flag;
  int ngood,nbad,good,bad,ipt;
  double value,param;
  int p[3],pflag[3];
  double *x,*xp[3];

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  for (int iface = 0; iface < 6; iface++) {
    dim = iface / 2;
    side = iface % 2;
    if (side == 0) value = boxlo[dim];
    else value = boxhi[dim];

    // check if tri straddles clipping plane
    // straddle = at least 2 pts on different sides (not on) clipping plane
    // case A: 1 pt on good side of clipping plane, 1 bad, 1 bad or on plane
    //   add 1 or 2 pts as needed on plane, keep as one tri
    // case B: 2 pts on good side of clipping plane, 1 on bad
    //   add 2 new points on plane, tri becomes trapezoid, so split into 2 tris
    // reallocate pts and tris if needed and reset xp
    // edge stores indices of pts already added, so don't create duplicate pts

    nedge = maxedge = 0;
    edge = NULL;
    
    m = ntri_old;
    n = ntri_new;
    for (i = 0; i < n; i++) {
      p[0] = tris[m].p1;
      p[1] = tris[m].p2;
      p[2] = tris[m].p3;
      xp[0] = pts[p[0]].x;
      xp[1] = pts[p[1]].x;
      xp[2] = pts[p[2]].x;

      if (side == 0) {
	if (xp[0][dim] < value) pflag[0] = BAD;
	else if (xp[0][dim] > value) pflag[0] = GOOD;
	else pflag[0] = NEITHER;
	if (xp[1][dim] < value) pflag[1] = BAD;
	else if (xp[1][dim] > value) pflag[1] = GOOD;
	else pflag[1] = NEITHER;
	if (xp[2][dim] < value) pflag[2] = BAD;
	else if (xp[2][dim] > value) pflag[2] = GOOD;
	else pflag[2] = NEITHER;
      } else {
	if (xp[0][dim] > value) pflag[0] = BAD;
	else if (xp[0][dim] < value) pflag[0] = GOOD;
	else pflag[0] = NEITHER;
	if (xp[1][dim] > value) pflag[1] = BAD;
	else if (xp[1][dim] < value) pflag[1] = GOOD;
	else pflag[1] = NEITHER;
	if (xp[2][dim] > value) pflag[2] = BAD;
	else if (xp[2][dim] < value) pflag[2] = GOOD;
	else pflag[2] = NEITHER;
      }

      nbad = ngood = 0;
      if (pflag[0] == BAD) nbad++;
      if (pflag[1] == BAD) nbad++;
      if (pflag[2] == BAD) nbad++;
      if (pflag[0] == GOOD) ngood++;
      if (pflag[1] == GOOD) ngood++;
      if (pflag[2] == GOOD) ngood++;

      if (nbad && ngood) {
	if (npoint_new+2 >= maxpoint) {
	  maxpoint += DELTA;
	  pts = (Surf::Point *) 
	    memory->srealloc(pts,maxpoint*sizeof(Surf::Point),"surf:pts");
	  xp[0] = pts[p[0]].x;
	  xp[1] = pts[p[1]].x;
	  xp[2] = pts[p[2]].x;
	}
	if (ntri_new == maxtri) {
	  maxtri += DELTA;
	  tris = (Surf::Tri *) 
	    memory->srealloc(tris,maxtri*sizeof(Surf::Tri),"surf:tris");
	}

	if (ngood == 1) {
	  if (pflag[0] == GOOD) good = 0;
	  else if (pflag[1] == GOOD) good = 1;
	  else if (pflag[2] == GOOD) good = 2;
	  if (pflag[0] == BAD) bad = 0;
	  else if (pflag[1] == BAD) bad = 1;
	  else if (pflag[2] == BAD) bad = 2;

	  ipt = find_edge(p[good],p[bad]);
	  if (ipt < 0) {
	    param = (value-xp[good][dim]) / (xp[bad][dim]-xp[good][dim]);
	    pts[npoint_new].x[0] = xp[good][0] + param*(xp[bad][0]-xp[good][0]);
	    pts[npoint_new].x[1] = xp[good][1] + param*(xp[bad][1]-xp[good][1]);
	    pts[npoint_new].x[2] = xp[good][2] + param*(xp[bad][2]-xp[good][2]);
	    pts[npoint_new].x[dim] = value;
	    ipt = npoint_new++;
	    add_edge(p[good],p[bad],ipt);
	  }

	  if (bad == 0) tris[m].p1 = ipt;
	  else if (bad == 1) tris[m].p2 = ipt;
	  else if (bad == 2) tris[m].p3 = ipt;

	  if (nbad == 2) {
	    if (bad == 0 && pflag[1] == BAD) bad = 1;
	    else bad = 2;

	    ipt = find_edge(p[good],p[bad]);
	    if (ipt < 0) {
	      param = (value-xp[good][dim]) / (xp[bad][dim]-xp[good][dim]);
	      pts[npoint_new].x[0] = 
		xp[good][0] + param*(xp[bad][0]-xp[good][0]);
	      pts[npoint_new].x[1] = 
		xp[good][1] + param*(xp[bad][1]-xp[good][1]);
	      pts[npoint_new].x[2] = 
		xp[good][2] + param*(xp[bad][2]-xp[good][2]);
	      pts[npoint_new].x[dim] = value;
	      ipt = npoint_new++;
	      add_edge(p[good],p[bad],ipt);
	    }

	    if (bad == 0) tris[m].p1 = ipt;
	    else if (bad == 1) tris[m].p2 = ipt;
	    else if (bad == 2) tris[m].p3 = ipt;
	  }

	} else if (ngood == 2) {
	  if (pflag[0] == GOOD) good = 0;
	  else if (pflag[1] == GOOD) good = 1;
	  else if (pflag[2] == GOOD) good = 2;
	  if (pflag[0] == BAD) bad = 0;
	  else if (pflag[1] == BAD) bad = 1;
	  else if (pflag[2] == BAD) bad = 2;

	  ipt = find_edge(p[good],p[bad]);
	  if (ipt < 0) {
	    param = (value-xp[good][dim]) / (xp[bad][dim]-xp[good][dim]);
	    pts[npoint_new].x[0] = xp[good][0] + param*(xp[bad][0]-xp[good][0]);
	    pts[npoint_new].x[1] = xp[good][1] + param*(xp[bad][1]-xp[good][1]);
	    pts[npoint_new].x[2] = xp[good][2] + param*(xp[bad][2]-xp[good][2]);
	    pts[npoint_new].x[dim] = value;
	    ipt = npoint_new++;
	    add_edge(p[good],p[bad],ipt);
	  }

	  if (bad == 0) tris[m].p1 = ipt;
	  else if (bad == 1) tris[m].p2 = ipt;
	  else if (bad == 2) tris[m].p3 = ipt;

	  memcpy(&tris[ntri_new],&tris[m],sizeof(Surf::Tri));
	  if (good == 0) tris[ntri_new].p1 = ipt;
	  else if (good == 1) tris[ntri_new].p2 = ipt;
	  else if (good == 2) tris[ntri_new].p3 = ipt;

	  if (good == 0 && pflag[1] == GOOD) good = 1;
	  else good = 2;

	  ipt = find_edge(p[good],p[bad]);
	  if (ipt < 0) {
	    param = (value-xp[good][dim]) / (xp[bad][dim]-xp[good][dim]);
	    pts[npoint_new].x[0] = xp[good][0] + param*(xp[bad][0]-xp[good][0]);
	    pts[npoint_new].x[1] = xp[good][1] + param*(xp[bad][1]-xp[good][1]);
	    pts[npoint_new].x[2] = xp[good][2] + param*(xp[bad][2]-xp[good][2]);
	    pts[npoint_new].x[dim] = value;
	    ipt = npoint_new++;
	    add_edge(p[good],p[bad],ipt);
	  }

	  if (bad == 0) tris[ntri_new].p1 = ipt;
	  else if (bad == 1) tris[ntri_new].p2 = ipt;
	  else if (bad == 2) tris[ntri_new].p3 = ipt;
	  ntri_new++;
	}
      }
      m++;
    }

    memory->destroy(edge);

    // project all points outside clipping plane to the plane

    m = npoint_old;
    for (i = 0; i < npoint_new; i++) {
      x = pts[m].x;
      if (side == 0 && x[dim] < value) x[dim] = value;
      if (side == 1 && x[dim] > value) x[dim] = value;
      m++;
    }

    // ptflag[I] = # of tris that include point I
    
    int *ptflag;
    memory->create(ptflag,npoint_new+npoint_old,"readsurf:ptflag");
    m = npoint_old;
    for (i = 0; i < npoint_new; i++) ptflag[m++] = 0;

    m = ntri_old;
    for (i = 0; i < ntri_new; i++) {
      ptflag[tris[m].p1]++;
      ptflag[tris[m].p2]++;
      ptflag[tris[m].p3]++;
      m++;
    }

    // triflag[I] = 1 if tri I can be removed b/c lies on clip plane
    // also decrement ptflag for pts in removed tri
    
    int *triflag;
    memory->create(triflag,ntri_new+ntri_old,"readsurf:triflag");
    m = ntri_old;
    for (i = 0; i < ntri_new; i++) triflag[m++] = 0;

    m = ntri_old;
    for (i = 0; i < ntri_new; i++) {
      if (pts[tris[m].p1].x[dim] == value && 
	  pts[tris[m].p2].x[dim] == value &&
	  pts[tris[m].p3].x[dim] == value) {
	triflag[m] = 1;
	ptflag[tris[m].p1]--;
	ptflag[tris[m].p2]--;
	ptflag[tris[m].p3]--;
      }
      m++;
    }

    // condense pts/tris to remove deleted points and tris
    // when delete points, set ptflag[I] = new index of point I for kept points
    // renumber point indices in tris using altered ptflag
    // reset npoint_new,ntri_new

    n = m = npoint_old;
    for (i = 0; i < npoint_new; i++) {
      if (ptflag[m]) {
	memcpy(&pts[n],&pts[m],sizeof(Surf::Point));
	ptflag[m] = n;
	n++;
      }
      m++;
    }
    npoint_new = n - npoint_old;
    
    n = m = ntri_old;
    for (i = 0; i < ntri_new; i++) {
      if (!triflag[m]) {
	memcpy(&tris[n],&tris[m],sizeof(Surf::Tri));
	tris[n].p1 = ptflag[tris[n].p1];
	tris[n].p2 = ptflag[tris[n].p2];
	tris[n].p3 = ptflag[tris[n].p3];
	n++;
      }
      m++;
    }
    ntri_new = n - ntri_old;

    // clean up

    memory->destroy(ptflag);
    memory->destroy(triflag);
  }

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  clipped to %d points\n",npoint_new);
      fprintf(screen,"  clipped to %d tris\n",ntri_new);
    }
    if (logfile) {
      fprintf(logfile,"  clipped to %d points\n",npoint_new);
      fprintf(logfile,"  clipped to %d tris\n",ntri_new);
    }
  }
}

/* ----------------------------------------------------------------------
   push all points to box boundary that are just inside
   delta = user-specified frac * (hi-lo)
   this avoids tiny clipped surf elements
------------------------------------------------------------------------- */

void ReadSurf::push_points_to_boundary(double frac)
{
  double *x;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  double xdelta = frac * (boxhi[0]-boxlo[0]);
  double ydelta = frac * (boxhi[1]-boxlo[1]);
  double zdelta = frac * (boxhi[2]-boxlo[2]);

  int n = npoint_old + npoint_new;
  for (int i = npoint_old; i < n; i++) {
    x = pts[i].x;
    if (x[0] >= boxlo[0] && x[0] <= boxhi[0]) {
      if (x[0]-boxlo[0] < xdelta) x[0] = boxlo[0];
      else if (boxhi[0]-x[0] < xdelta) x[0] = boxhi[0];
    }
    if (x[1] >= boxlo[1] && x[1] <= boxhi[1]) {
      if (x[1]-boxlo[1] < ydelta) x[1] = boxlo[1];
      else if (boxhi[1]-x[1] < ydelta) x[1] = boxhi[1];
    }
    if (dimension == 2) continue;
    if (x[2] >= boxlo[2] && x[2] <= boxhi[2]) {
      if (x[2]-boxlo[2] < zdelta) x[2] = boxlo[2];
      else if (boxhi[2]-x[2] < zdelta) x[2] = boxhi[2];
    }
  }
}

/* ----------------------------------------------------------------------
   check if all new points are inside or on surface of global simulation box
------------------------------------------------------------------------- */

void ReadSurf::check_point_inside()
{
  double *x;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  int m = npoint_old;
  int nbad = 0;
  for (int i = 0; i < npoint_new; i++) {
    x = pts[m].x;
    if (x[0] < boxlo[0] || x[0] > boxhi[0] ||
	x[1] < boxlo[1] || x[1] > boxhi[1] ||
	x[2] < boxlo[2] || x[2] > boxhi[2]) nbad++;
    m++;
  }

  if (nbad) {
    char str[128];
    sprintf(str,"%d read_surf points are not inside simulation box",
	    nbad);
    error->all(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   check that new points are each end point of exactly 2 new lines
   exception: not required of point on simulation box surface
------------------------------------------------------------------------- */

void ReadSurf::check_watertight_2d()
{
  int p1,p2;

  // count[I] = # of lines that vertex I is part of

  int *count;
  memory->create(count,npoint_new,"readsurf:count");
  for (int i = 0; i < npoint_new; i++) count[i] = 0;

  int ndup = 0;
  int m = nline_old;
  for (int i = 0; i < nline_new; i++) {
    p1 = lines[m].p1 - npoint_old;
    p2 = lines[m].p2 - npoint_old;
    count[p1]++;
    count[p2]++;
    if (count[p1] > 2) ndup++;
    if (count[p2] > 2) ndup++;
    m++;
  }
  
  if (ndup) {
    char str[128];
    sprintf(str,"Surface check failed with %d duplicate points",ndup);
    error->all(FLERR,str);
  }

  // check that all counts are 2
  // allow for exception if point on box surface

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  int nbad = 0;
  for (int i = 0; i < npoint_new; i++) {
    if (count[i] == 0) nbad++;
    else if (count[i] == 1) {
      if (!Geometry::point_on_hex(pts[i+npoint_old].x,boxlo,boxhi)) nbad++;
    }
  }

  if (nbad) {
    char str[128];
    sprintf(str,"Surface check failed with %d unmatched points",nbad);
    error->all(FLERR,str);
  }

  // clean up

  memory->destroy(count);
}

/* ----------------------------------------------------------------------
   check new directed triangle edges
   must be unique and match exactly one inverted edge
   exception: not required of triangle edge on simulation box surface
------------------------------------------------------------------------- */

void ReadSurf::check_watertight_3d()
{
  // hash directed edges of all triangles
  // key = directed edge, value = # of times it appears in any triangle
  // NOTE: could prealloc hash to correct size here

#ifdef SPARTA_MAP
  std::map<bigint,int> hash;
  std::map<bigint,int>::iterator it;
#elif SPARTA_UNORDERED_MAP
  std::unordered_map<bigint,int> hash;
  std::unordered_map<bigint,int>::iterator it;
#else
  std::tr1::unordered_map<bigint,int> hash;
  std::tr1::unordered_map<bigint,int>::iterator it;
#endif

  // insert each edge into hash
  // should appear once in each direction
  // error if any duplicate edges

  bigint p1,p2,p3,key;

  int ndup = 0;
  int m = ntri_old;
  for (int i = 0; i < ntri_new; i++) {
    p1 = tris[m].p1;
    p2 = tris[m].p2;
    p3 = tris[m].p3;
    key = (p1 << 32) | p2;
    if (hash.find(key) != hash.end()) ndup++;
    else hash[key] = 0;
    key = (p2 << 32) | p3;
    if (hash.find(key) != hash.end()) ndup++;
    else hash[key] = 0;
    key = (p3 << 32) | p1;
    if (hash.find(key) != hash.end()) ndup++;
    else hash[key] = 0;
    m++;
  }

  if (ndup) {
    char str[128];
    sprintf(str,"Surface check failed with %d duplicate edges",ndup);
    error->all(FLERR,str);
  }

  // check that each edge has an inverted match
  // allow for exception if edge on box surface

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  int nbad = 0;
  for (it = hash.begin(); it != hash.end(); ++it) {
    p1 = it->first >> 32;
    p2 = it->first & MAXSMALLINT;
    key = (p2 << 32) | p1;
    if (hash.find(key) == hash.end()) {
      if (!Geometry::point_on_hex(pts[p1].x,boxlo,boxhi) ||
          !Geometry::point_on_hex(pts[p2].x,boxlo,boxhi)) nbad++;
    }
  }

  if (nbad) {
    char str[128];
    sprintf(str,"Surface check failed with %d unmatched edges",nbad);
    error->all(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   check norms of new adjacent lines
   error if dot product of 2 norms is -1 -> infinitely thin surf
   warn if closer than EPSILON_NORM to -1
------------------------------------------------------------------------- */

void ReadSurf::check_neighbor_norm_2d()
{
  int p1,p2;

  // count[I] = # of lines that vertex I is part of

  int *count;
  int **p2e;
  memory->create(count,npoint_new,"readsurf:count");
  memory->create(p2e,npoint_new,2,"readsurf:count");
  for (int i = 0; i < npoint_new; i++) count[i] = 0;

  int m = nline_old;
  for (int i = 0; i < nline_new; i++) {
    p1 = lines[m].p1 - npoint_old;
    p2 = lines[m].p2 - npoint_old;
    p2e[p1][count[p1]++] = i;
    p2e[p2][count[p2]++] = i;
    m++;
  }
  
  // check that norms of adjacent lines are not in opposite directions

  double dot;
  double *norm1,*norm2;

  int nerror = 0;
  int nwarn = 0;
  for (int i = 0; i < npoint_new; i++) {
    if (count[i] == 1) continue;
    norm1 = lines[p2e[i][0]].norm;
    norm2 = lines[p2e[i][1]].norm;
    dot = MathExtra::dot3(norm1,norm2);
    if (dot <= -1.0) nerror++;
    else if (dot < -1.0+EPSILON_NORM) nwarn++;

  }

  if (nerror) {
    char str[128];
    sprintf(str,"Surface check failed with %d "
            "infinitely thin line pairs",nerror);
    error->all(FLERR,str);
  }
  if (nwarn) {
    char str[128];
    sprintf(str,"Surface check found %d "
            "nearly infinitely thin line pairs",nwarn);
    if (me == 0) error->warning(FLERR,str);
  }

  // clean up

  memory->destroy(count);
  memory->destroy(p2e);
}

/* ----------------------------------------------------------------------
   check norms of new adjacent triangles
   error if dot product of 2 norms is -1 -> infinitely thin surf
   warn if closer than EPSILON_NORM to -1
------------------------------------------------------------------------- */

void ReadSurf::check_neighbor_norm_3d()
{
  // hash directed edges of all triangles
  // key = directed edge, value = triangle it is part of
  // NOTE: could prealloc hash to correct size here

#ifdef SPARTA_MAP
  std::map<bigint,int> hash;
  std::map<bigint,int>::iterator it;
#elif SPARTA_UNORDERED_MAP
  std::unordered_map<bigint,int> hash;
  std::unordered_map<bigint,int>::iterator it;
#else
  std::tr1::unordered_map<bigint,int> hash;
  std::tr1::unordered_map<bigint,int>::iterator it;
#endif

  // insert each edge into hash with triangle as value

  bigint p1,p2,p3,key;

  int m = ntri_old;
  for (int i = 0; i < ntri_new; i++) {
    p1 = tris[m].p1;
    p2 = tris[m].p2;
    p3 = tris[m].p3;
    key = (p1 << 32) | p2;
    hash[key] = m;
    key = (p2 << 32) | p3;
    hash[key] = m;
    key = (p3 << 32) | p1;
    hash[key] = m;
    m++;
  }

  // check that norms of adjacent triangles are not in opposite directions

  double dot;
  double *norm1,*norm2;

  int nerror = 0;
  int nwarn = 0;
  for (it = hash.begin(); it != hash.end(); ++it) {
    p1 = it->first >> 32;
    p2 = it->first & MAXSMALLINT;
    key = (p2 << 32) | p1;
    if (hash.find(key) == hash.end()) continue;
    norm1 = tris[it->second].norm;
    norm2 = tris[hash[key]].norm;
    dot = MathExtra::dot3(norm1,norm2);
    if (dot <= -1.0) nerror++;
    else if (dot < -1.0+EPSILON_NORM) nwarn++;
  }

  if (nerror) {
    char str[128];
    sprintf(str,"Surface check failed with %d "
            "infinitely thin triangle pairs",nerror);
    error->all(FLERR,str);
  }
  if (nwarn) {
    char str[128];
    sprintf(str,"Surface check found %d "
            "nearly infinitely thin triangle pairs",nwarn);
    if (me == 0) error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   check nearness of all points to other lines in same cell
   error if point is on line, including duplicate point
   warn if closer than EPSILON_GRID = fraction of grid cell size
   NOTE: this can miss a close point/line pair in 2 different grid cells
------------------------------------------------------------------------- */

void ReadSurf::check_point_near_surf_2d()
{
  int i,j,n,p1,p2;
  int *csurfs;
  double side,epssq;
  double delta[3];
  double *lo,*hi;
  Surf::Line *line;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  int nerror = 0;
  int nwarn = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    n = cells[icell].nsurf;
    if (n == 0) continue;

    lo = cells[icell].lo;
    hi = cells[icell].hi;
    side = MIN(hi[0]-lo[0],hi[1]-lo[1]);
    epssq = (EPSILON_GRID*side) * (EPSILON_GRID*side);

    csurfs = cells[icell].csurfs;
    for (i = 0; i < n; i++) {
      line = &lines[csurfs[i]];
      for (j = 0; j < n; j++) {
        if (i == j) continue;
        p1 = lines[csurfs[j]].p1;
        p2 = lines[csurfs[j]].p2;
        point_line_compare(p1,line,epssq,nerror,nwarn);
        point_line_compare(p2,line,epssq,nerror,nwarn);
      }
    }
  }

  int all;
  MPI_Allreduce(&nerror,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,"Surface check failed with %d points on lines",all);
    error->all(FLERR,str);
  }

  MPI_Allreduce(&nwarn,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,"Surface check found %d points nearly on lines",all);
    if (me == 0) error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   check nearness of all points to other triangles in same cell
   error if point is on triangle, including duplicate point
   warn if closer than EPSILON_GRID = fraction of grid cell size
   NOTE: this can miss a close point/triangle pair in 2 different grid cells
------------------------------------------------------------------------- */

void ReadSurf::check_point_near_surf_3d()
{
  int i,j,n,p1,p2,p3;
  int *csurfs;
  double side,epssq;
  double delta[3];
  double *lo,*hi;
  Surf::Tri *tri;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  int nerror = 0;
  int nwarn = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    n = cells[icell].nsurf;
    if (n == 0) continue;

    lo = cells[icell].lo;
    hi = cells[icell].hi;
    side = MIN(hi[0]-lo[0],hi[1]-lo[1]);
    side = MIN(side,hi[2]-lo[2]);
    epssq = (EPSILON_GRID*side) * (EPSILON_GRID*side);

    csurfs = cells[icell].csurfs;
    for (i = 0; i < n; i++) {
      tri = &tris[csurfs[i]];
      for (j = 0; j < n; j++) {
        if (i == j) continue;
        p1 = tris[csurfs[j]].p1;
        p2 = tris[csurfs[j]].p2;
        p3 = tris[csurfs[j]].p3;
        point_tri_compare(p1,tri,epssq,nerror,nwarn,icell,csurfs[i],csurfs[j]);
        point_tri_compare(p2,tri,epssq,nerror,nwarn,icell,csurfs[i],csurfs[j]);
        point_tri_compare(p3,tri,epssq,nerror,nwarn,icell,csurfs[i],csurfs[j]);
      }
    }
  }

  int all;
  MPI_Allreduce(&nerror,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,"Surface check failed with %d points on triangles",all);
    error->all(FLERR,str);
  }

  MPI_Allreduce(&nwarn,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,"Surface check found %d points nearly on triangles",all);
    if (me == 0) error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   compute distance bewteen a point and line
   just return if point is an endpoint of line
   increment nerror if point on line
   increment nwarn if point is within epssq distance of line
------------------------------------------------------------------------- */

void ReadSurf::point_line_compare(int i, Surf::Line *line, double epssq, 
                                  int &nerror, int &nwarn)
{
  if (i == line->p1 || i == line->p2) return;
  double rsq = 
    Geometry::distsq_point_line(pts[i].x,pts[line->p1].x,pts[line->p2].x);
  if (rsq == 0.0) nerror++;
  else if (rsq < epssq) nwarn++;
}

/* ----------------------------------------------------------------------
   compute distance bewteen a point and triangle
   just return if point is an endpoint of triangle
   increment nerror if point on triangle
   increment nwarn if point is within epssq distance of triangle
------------------------------------------------------------------------- */

void ReadSurf::point_tri_compare(int i, Surf::Tri *tri, double epssq, 
                                 int &nerror, int &nwarn, 
                                 int icell, int ci, int cj)
{
  if (i == tri->p1 || i == tri->p2 || i == tri->p3) return;
  double rsq = 
    Geometry::distsq_point_tri(pts[i].x,
                               pts[tri->p1].x,pts[tri->p2].x,pts[tri->p3].x,
                               tri->norm);
  if (rsq == 0.0) nerror++;
  else if (rsq < epssq) nwarn++;
}

/* ----------------------------------------------------------------------
   find edge (I,J) in either order within edge list
   return index of new point already created along the edge
------------------------------------------------------------------------- */

int ReadSurf::find_edge(int i, int j)
{
  for (int m = 0; m < nedge; m++) {
    if (i == edge[m][0] && j == edge[m][1]) return edge[m][2];
    if (i == edge[m][1] && j == edge[m][0]) return edge[m][2];
  }
  return -1;
}

/* ----------------------------------------------------------------------
   add edge (I,J) to edge list with index of new pt M created along that edge
   grow edge list if necessary
------------------------------------------------------------------------- */

void ReadSurf::add_edge(int i, int j, int m)
{
  if (nedge == maxedge) {
    maxedge += DELTA;
    memory->grow(edge,maxedge,3,"readsurf:edge");
  }

  edge[nedge][0] = i;
  edge[nedge][1] = j;
  edge[nedge][2] = m;
  nedge++;
}

/* ----------------------------------------------------------------------
   return shortest line length
------------------------------------------------------------------------- */

double ReadSurf::shortest_line()
{
  double len = BIG;
  int m = nline_old;
  for (int i = 0; i < nline_new; i++) {
    len = MIN(len,surf->line_size(m));
    m++;
  }
  return len;
}

/* ----------------------------------------------------------------------
   return shortest tri edge and smallest tri area
------------------------------------------------------------------------- */

void ReadSurf::smallest_tri(double &len, double &area)
{
  double lenone,areaone;

  len = area = BIG;
  int m = ntri_old;
  for (int i = 0; i < ntri_new; i++) {
    areaone = surf->tri_size(m,lenone);
    len = MIN(len,lenone);
    area = MIN(area,areaone);
    m++;
  }
}

/* ----------------------------------------------------------------------
   proc 0 opens data file
   test if gzipped
------------------------------------------------------------------------- */

void ReadSurf::open(char *file)
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
   grab next keyword
   read lines until one is non-blank
   keyword is all text on line w/out leading & trailing white space
   read one additional line (assumed blank)
   if any read hits EOF, set keyword to empty
   if first = 1, line variable holds non-blank line that ended header
------------------------------------------------------------------------- */

void ReadSurf::parse_keyword(int first)
{
  int eof = 0;

  // proc 0 reads upto non-blank line plus 1 following line
  // eof is set to 1 if any read hits end-of-file

  if (me == 0) {
    if (!first) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
    }
    while (eof == 0 && strspn(line," \t\n\r") == strlen(line)) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
    }
    if (fgets(buffer,MAXLINE,fp) == NULL) eof = 1;
  }

  // if eof, set keyword empty and return

  MPI_Bcast(&eof,1,MPI_INT,0,world);
  if (eof) {
    keyword[0] = '\0';
    return;
  }

  // bcast keyword line to all procs

  int n;
  if (me == 0) n = strlen(line) + 1;
  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);

  // copy non-whitespace portion of line into keyword

  int start = strspn(line," \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t' 
	 || line[stop] == '\n' || line[stop] == '\r') stop--;
  line[stop+1] = '\0';
  strcpy(keyword,&line[start]);
}

/* ----------------------------------------------------------------------
   count and return words in a single line
   make copy of line before using strtok so as not to change line
   trim anything from '#' onward
------------------------------------------------------------------------- */

int ReadSurf::count_words(char *line)
{
  int n = strlen(line) + 1;
  char *copy = (char *) memory->smalloc(n*sizeof(char),"copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->sfree(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->sfree(copy);
  return n;
}

/* ----------------------------------------------------------------------
   unneeded code for now
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   check if any pair of points is closer than epsilon
   done in O(N) manner by binning twice with offset bins
  //NOTE: check if bins allow for surf points
  //NOTE: or maybe should ignore surf points in this test, since clip
  //      could put them very close together
------------------------------------------------------------------------- */

/* 

void ReadSurf::check_point_pairs()
{
  int i,j,k,m,n,ix,iy,iz;
  double dx,dy,dz,rsq;
  double origin[3];

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  // epsilon = EPSILON fraction of shortest box length
  // epssq = epsilon squared

  double epsilon = MIN(domain->xprd,domain->yprd);
  if (dimension == 3) epsilon = MIN(epsilon,domain->zprd);
  epsilon *= EPSILON;
  double epssq = epsilon * epsilon;

  // goal: N roughly cubic bins where N = # of new points
  // nbinxyz = # of bins in each dimension
  // xyzbin = bin size in each dimension
  // for 2d, nbinz = 1
  // after setting bin size, add 1 to nbinxyz
  // this allows for 2nd binning via offset origin

  int nbinx,nbiny,nbinz;
  double xbin,ybin,zbin;
  double xbininv,ybininv,zbininv;
  
  if (dimension == 2) {
    double vol_per_point = domain->xprd * domain->yprd / npoint_new;
    xbin = ybin = sqrt(vol_per_point);
    nbinx = static_cast<int> (domain->xprd / xbin);
    nbiny = static_cast<int> (domain->yprd / ybin);
    if (nbinx == 0) nbinx = 1;
    if (nbiny == 0) nbiny = 1;
    nbinz = 1;
  } else {
    double vol_per_point = domain->xprd * domain->yprd * domain->zprd / 
      npoint_new;
    xbin = ybin = zbin = pow(vol_per_point,1.0/3.0);
    nbinx = static_cast<int> (domain->xprd / xbin);
    nbiny = static_cast<int> (domain->yprd / ybin);
    nbinz = static_cast<int> (domain->zprd / zbin);
    if (nbinx == 0) nbinx = 1;
    if (nbiny == 0) nbiny = 1;
    if (nbinz == 0) nbinz = 1;
  }

  xbin = domain->xprd / nbinx;
  ybin = domain->yprd / nbiny;
  zbin = domain->zprd / nbinz;
  xbininv = 1.0/xbin;
  ybininv = 1.0/ybin;
  zbininv = 1.0/zbin;

  if (nbinx > 1) nbinx++;
  if (nbiny > 1) nbiny++;
  if (nbinz > 1) nbinz++;

  // binhead[I][J][K] = point index of 1st point in bin IJK, -1 if none
  // bin[I] = index of next point in same bin as point I, -1 if last

  int ***binhead;
  memory->create(binhead,nbinx,nbiny,nbinz,"readsurf:binhead");
  int *bin;
  memory->create(bin,npoint_new,"readsurf:bin");

  // 1st binning = bins aligned with global box boundaries

  for (i = 0; i < nbinx; i++)
    for (j = 0; j < nbiny; j++)
      for (k = 0; k < nbinz; k++)
	binhead[i][j][k] = -1;

  origin[0] = boxlo[0];
  origin[1] = boxlo[1];
  origin[2] = boxlo[2];

  m = npoint_old;
  for (i = 0; i < npoint_new; i++) {
    ix = static_cast<int> ((pts[m].x[0] - origin[0]) * xbininv);
    iy = static_cast<int> ((pts[m].x[1] - origin[1]) * ybininv);
    iz = static_cast<int> ((pts[m].x[2] - origin[2]) * zbininv);
    bin[m] = binhead[ix][iy][iz];
    binhead[ix][iy][iz] = m;
    m++;
  }

  // check distances for all pairs of particles in same bin

  int nbad = 0;

  for (i = 0; i < nbinx; i++)
    for (j = 0; j < nbiny; j++)
      for (k = 0; k < nbinz; k++) {
	m = binhead[i][j][k];
	while (m >= 0) {
	  n = bin[m];
	  while (n >= 0) {
	    dx = pts[m].x[0] - pts[n].x[0];
	    dy = pts[m].x[1] - pts[n].x[1];
	    dz = pts[m].x[2] - pts[n].x[2];
	    rsq = dx*dx + dy*dy + dz*dz;
	    if (rsq < epssq) nbad++;
	    n = bin[n];
	  }
	  m = bin[m];
	}
      }

  if (nbad) {
    char str[128];
    sprintf(str,"%d read_surf point pairs are too close",nbad);
    error->all(FLERR,str);
  }

  // 2nd binning = bins offset by 1/2 binsize wrt global box boundaries
  // do not offset bin origin in a dimension with only 1 bin

  for (i = 0; i < nbinx; i++)
    for (j = 0; j < nbiny; j++)
      for (k = 0; k < nbinz; k++)
	binhead[i][j][k] = -1;

  origin[0] = boxlo[0] - 0.5*xbin;
  origin[1] = boxlo[1] - 0.5*ybin;
  origin[2] = boxlo[2] - 0.5*zbin;
  if (nbinx == 1) origin[0] = boxlo[0];
  if (nbiny == 1) origin[1] = boxlo[1];
  if (nbinz == 1) origin[2] = boxlo[2];

  m = npoint_old;
  for (i = 0; i < npoint_new; i++) {
    ix = static_cast<int> ((pts[m].x[0] - origin[0]) * xbininv);
    iy = static_cast<int> ((pts[m].x[1] - origin[1]) * ybininv);
    iz = static_cast<int> ((pts[m].x[2] - origin[2]) * zbininv);
    bin[m] = binhead[ix][iy][iz];
    binhead[ix][iy][iz] = m;
    m++;
  }

  // check distances for all pairs of particles in same bin

  nbad = 0;

  for (i = 0; i < nbinx; i++)
    for (j = 0; j < nbiny; j++)
      for (k = 0; k < nbinz; k++) {
	m = binhead[i][j][k];
	while (m >= 0) {
	  n = bin[m];
	  while (n >= 0) {
	    dx = pts[m].x[0] - pts[n].x[0];
	    dy = pts[m].x[1] - pts[n].x[1];
	    dz = pts[m].x[2] - pts[n].x[2];
	    rsq = dx*dx + dy*dy + dz*dz;
	    if (rsq < epssq) nbad++;
	    n = bin[n];
	  }
	  m = bin[m];
	}
      }

  if (nbad) {
    char str[128];
    sprintf(str,"%d read_surf point pairs are too close",nbad);
    error->all(FLERR,str);
  }

  // clean up

  memory->destroy(binhead);
  memory->destroy(bin);
}

*/
