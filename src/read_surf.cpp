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
#include "write_surf.h"
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

  surf->exist = 1;
  dim = domain->dimension;

  if (narg < 1) error->all(FLERR,"Illegal read_surf command");

  // read header info

  if (me == 0) {
    if (screen) fprintf(screen,"Reading surf file ...\n");
    open(arg[0]);
  }

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  header();

  // create pts,lines,tris data structures

  pts = (Point *) memory->smalloc(npoint*sizeof(Point),"readsurf:pts");
  lines = (Line *) memory->smalloc(nline*sizeof(Line),"readsurf:lines");
  tris = (Tri *) memory->smalloc(ntri*sizeof(Tri),"readsurf:tris");

  maxpoint = npoint;
  maxline = nline;
  maxtri = ntri;

  // read and store Points and Lines/Tris sections

  parse_keyword(1);
  if (strcmp(keyword,"Points") != 0)
    error->all(FLERR,
	       "Read_surf did not find points section of surf file");
  read_points();

  parse_keyword(0);
  if (dim == 2) {
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
  // store optional keyword for file output

  origin[0] = origin[1] = origin[2] = 0.0;
  int grouparg = 0;
  int typeadd = 0;
  int partflag = NONE;
  int filearg = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"origin") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Invalid read_surf command");
      double ox = atof(arg[iarg+1]);
      double oy = atof(arg[iarg+2]);
      double oz = atof(arg[iarg+3]);
      if (dim == 2 && oz != 0.0) 
	error->all(FLERR,"Invalid read_surf geometry transformation "
		   "for 2d simulation");
      origin[0] = ox;
      origin[1] = oy;
      origin[2] = oz;
      iarg += 4;
    } else if (strcmp(arg[iarg],"trans") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Invalid read_surf command");
      double dx = input->numeric(FLERR,arg[iarg+1]);
      double dy = input->numeric(FLERR,arg[iarg+2]);
      double dz = input->numeric(FLERR,arg[iarg+3]);
      if (dim == 2 && dz != 0.0) 
	error->all(FLERR,"Invalid read_surf geometry transformation "
		   "for 2d simulation");
      origin[0] += dx;
      origin[1] += dy;
      origin[2] += dz;
      translate(dx,dy,dz);
      iarg += 4;
    } else if (strcmp(arg[iarg],"atrans") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Invalid read_surf command");
      double ax = input->numeric(FLERR,arg[iarg+1]);
      double ay = input->numeric(FLERR,arg[iarg+2]);
      double az = input->numeric(FLERR,arg[iarg+3]);
      if (dim == 2 && az != 0.0) 
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
      double fx = input->numeric(FLERR,arg[iarg+1]);
      double fy = input->numeric(FLERR,arg[iarg+2]);
      double fz = input->numeric(FLERR,arg[iarg+3]);
      if (dim == 2 && fz != 0.5) 
	error->all(FLERR,"Invalid read_surf geometry transformation "
		   "for 2d simulation");
      double ax = domain->boxlo[0] + fx*domain->xprd;
      double ay = domain->boxlo[1] + fy*domain->yprd;
      double az;
      if (dim == 3) az = domain->boxlo[2] + fz*domain->zprd;
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
      double sx = input->numeric(FLERR,arg[iarg+1]);
      double sy = input->numeric(FLERR,arg[iarg+2]);
      double sz = input->numeric(FLERR,arg[iarg+3]);
      if (dim == 2 && sz != 1.0) 
	error->all(FLERR,"Invalid read_surf geometry transformation "
		   "for 2d simulation");
      scale(sx,sy,sz);
      iarg += 4;
    } else if (strcmp(arg[iarg],"rotate") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Invalid read_surf command");
      double theta = input->numeric(FLERR,arg[iarg+1]);
      double rx = input->numeric(FLERR,arg[iarg+2]);
      double ry = input->numeric(FLERR,arg[iarg+3]);
      double rz = input->numeric(FLERR,arg[iarg+4]);
      if (dim == 2 && (rx != 0.0 || ry != 0.0 || rz != 1.0))
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
          frac = input->numeric(FLERR,arg[iarg+1]);
          if (frac < 0.0 || frac >= 0.5) 
            error->all(FLERR,"Invalid read_surf command");
          iarg++;
        }
      }
      if (frac > 0.0) push_points_to_boundary(frac);
      if (dim == 2) clip2d();
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

    } else if (strcmp(arg[iarg],"particle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid read_surf command");
      if (strcmp(arg[iarg+1],"none") == 0) partflag = NONE;
      else if (strcmp(arg[iarg+1],"check") == 0) partflag = CHECK;
      else if (strcmp(arg[iarg+1],"keep") == 0) partflag = KEEP;
      else error->all(FLERR,"Invalid read_surf command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid read_surf command");
      filearg = iarg+1;
      iarg += 2;

    } else error->all(FLERR,"Invalid read_surf command");
  }

  // error test on particles

  if (particle->exist && partflag == NONE)
    error->all(FLERR,"Using read_surf particle none when particles exist");

  // if specified, apply group and typeadd keywords
  // these reset per-element mask/type info

  if (grouparg) {
    int igroup = surf->find_group(arg[grouparg]);
    if (igroup < 0) igroup = surf->add_group(arg[grouparg]);
    int groupbit = surf->bitmask[igroup];
    if (dim == 2) {
      for (int i = 0; i < nline; i++) lines[i].mask |= groupbit;
    } else {
      for (int i = 0; i < ntri; i++) tris[i].mask |= groupbit;
    }
  }

  if (typeadd) {
    if (dim == 2) {
      for (int i = 0; i < nline; i++) lines[i].type += typeadd;
    } else {
      for (int i = 0; i < ntri; i++) tris[i].type += typeadd;
    }
  }

  // extent of surfs after geometric transformations
  // compute sizes of smallest surface elements

  double extent[3][2];
  extent[0][0] = extent[1][0] = extent[2][0] = BIG;
  extent[0][1] = extent[1][1] = extent[2][1] = -BIG;

  for (int i = 0; i < npoint; i++) {
    extent[0][0] = MIN(extent[0][0],pts[i].x[0]);
    extent[0][1] = MAX(extent[0][1],pts[i].x[0]);
    extent[1][0] = MIN(extent[1][0],pts[i].x[1]);
    extent[1][1] = MAX(extent[1][1],pts[i].x[1]);
    extent[2][0] = MIN(extent[2][0],pts[i].x[2]);
    extent[2][1] = MAX(extent[2][1],pts[i].x[2]);
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

  // write out new surf file if requested
  // do this before assigning surfs to grid cells, in case an error occurs

  if (filearg) {
    WriteSurf *wf = new WriteSurf(sparta);
    if (comm->me == 0) {
      FILE *fp = fopen(arg[filearg],"w");
      if (!fp) {
	char str[128];
	sprintf(str,"Cannot open surface file %s",arg[0]);
	error->one(FLERR,str);
      }
      wf->write_file(fp);
      fclose(fp);
    }
    delete wf;
  }

  // add read-in pts/lines/tris to Surf data structures

  if (dim == 2) {
    nline_old = surf->nline;
    nline_new = nline_old + nline;
    Surf::Line *newlines = surf->lines;
    newlines = (Surf::Line *) 
      memory->srealloc(newlines,nline_new*sizeof(Surf::Line),"surf:lines");

    int m = nline_old;
    for (int i = 0; i < nline; i++) {
      newlines[m].id = m+1;
      newlines[m].type = lines[i].type;
      newlines[m].mask = lines[i].mask;
      newlines[m].isc = newlines[m].isr = -1;
      memcpy(newlines[m].p1,pts[lines[i].p1].x,3*sizeof(double));
      memcpy(newlines[m].p2,pts[lines[i].p2].x,3*sizeof(double));
      m++;
    }

    surf->lines = newlines;

  } else if (dim == 3) {
    ntri_old = surf->ntri;
    ntri_new = ntri_old + ntri;
    Surf::Tri *newtris = surf->tris;
    newtris = (Surf::Tri *) 
      memory->srealloc(newtris,ntri_new*sizeof(Surf::Tri),"surf:tris");

    int m = ntri_old;
    for (int i = 0; i < ntri; i++) {
      newtris[m].id = m+1;
      newtris[m].type = tris[i].type;
      newtris[m].mask = tris[i].mask;
      newtris[m].isc = newtris[m].isr = -1;
      memcpy(newtris[m].p1,pts[tris[i].p1].x,3*sizeof(double));
      memcpy(newtris[m].p2,pts[tris[i].p2].x,3*sizeof(double));
      memcpy(newtris[m].p3,pts[tris[i].p3].x,3*sizeof(double));
      m++;
    }

    surf->tris = newtris;
  }

  // compute normals of new lines or triangles

  if (dim == 2) surf->compute_line_normal(nline_old,nline_new);
  else surf->compute_tri_normal(ntri_old,ntri_new);

  // error check on new points,lines,tris
  // all points must be inside or on surface of simulation box
  // NOTE: need to change this call

  //surf->check_point_inside(npoint_old,npoint_new);

  // -----------------------
  // map surfs to grid cells
  // -----------------------

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // sort particles

  if (particle->exist) particle->sort();

  // make list of surf elements I own
  // assign surfs to grid cells
  // error checks to flag bad surfs

  surf->setup_surf();

  grid->unset_neighbors();
  grid->remove_ghosts();

  if (particle->exist && grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
	grid->combine_split_cell_particles(icell,1);
  }

  grid->clear_surf();

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // error checks that can be done before surfs are mapped to grid cells

  if (dim == 2) {
    surf->check_watertight_2d(nline_old);
    check_neighbor_norm_2d();
  } else {
    surf->check_watertight_3d(ntri_old);
    check_neighbor_norm_3d();
  }

  // can now free local copy of read-in surfs
  // last use was in check_neighbor() methods

  memory->sfree(pts);
  memory->sfree(lines);
  memory->sfree(tris);

  MPI_Barrier(world);
  double time4 = MPI_Wtime();

  // map surfs to grid cells then error check
  // check done on per-grid-cell basis, too expensive to do globally

  grid->surf2grid(1);

  if (dim == 2) check_point_near_surf_2d();
  else check_point_near_surf_3d();

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

  // DEBUG
  //grid->debug();

  MPI_Barrier(world);
  double time7 = MPI_Wtime();

  // remove particles in any cell that is now INSIDE or has new surfs
  // reassign particles in split cells to sub cell owner
  // compress particles if any flagged for deletion

  bigint ndeleted;
  if (particle->exist) {
    Grid::ChildCell *cells = grid->cells;
    Grid::ChildInfo *cinfo = grid->cinfo;
    int nglocal = grid->nlocal;
    int delflag = 0;

    for (int icell = 0; icell < nglocal; icell++) {
      if (cinfo[icell].type == INSIDE) {
	if (partflag == KEEP) 
          error->one(FLERR,"Particles are inside new surfaces");
	if (cinfo[icell].count) delflag = 1;
	particle->remove_all_from_cell(cinfo[icell].first);
	cinfo[icell].count = 0;
	cinfo[icell].first = -1;
	continue;
      }
      if (cells[icell].nsurf && cells[icell].nsplit >= 1) {
	int nsurf = cells[icell].nsurf;
	int *csurfs = cells[icell].csurfs;
	int m;
	if (dim == 2) {
	  for (m = 0; m < nsurf; m++) {
	    if (csurfs[m] >= nline_old) break;
	  }
	} else {
	  for (m = 0; m < nsurf; m++) {
	    if (csurfs[m] >= ntri_old) break;
	  }
	}
	if (m < nsurf && partflag == CHECK) {
	  if (cinfo[icell].count) delflag = 1;
	  particle->remove_all_from_cell(cinfo[icell].first);
	  cinfo[icell].count = 0;
	  cinfo[icell].first = -1;
	}
      }
      if (cells[icell].nsplit > 1)
	grid->assign_split_cell_particles(icell);
    }
    int nlocal_old = particle->nlocal;
    if (delflag) particle->compress_rebalance();
    bigint delta = nlocal_old - particle->nlocal;
    MPI_Allreduce(&delta,&ndeleted,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  }

  MPI_Barrier(world);
  double time8 = MPI_Wtime();

  double time_total = time6-time1;

  if (comm->me == 0) {
    if (screen) {
      if (particle->exist)
	fprintf(screen,"  " BIGINT_FORMAT " deleted particles\n",ndeleted);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  read/sort/check/surf2grid/ghost/inout/particle percent = "
              "%g %g %g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total,100.0*(time7-time6)/time_total,
              100.0*(time8-time7)/time_total);
    }
    if (logfile) {
      if (particle->exist)
	fprintf(logfile,"  " BIGINT_FORMAT " deleted particles\n",ndeleted);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  read/sort/check/surf2grid/ghost/inout/particle percent = "
              "%g %g %g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total,100.0*(time7-time6)/time_total,
              100.0*(time8-time7)/time_total);
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

    // search line for header keyword and set corresponding variable

    if (strstr(line,"points")) sscanf(line,"%d",&npoint);
    else if (strstr(line,"lines")) {
      if (dim == 3) 
	error->all(FLERR,"Surf file cannot contain lines for 3d simulation");
      sscanf(line,"%d",&nline);
    } else if (strstr(line,"triangles")) {
      if (dim == 2) 
	error->all(FLERR,
		   "Surf file cannot contain triangles for 2d simulation");
      sscanf(line,"%d",&ntri);
    } else break;
  }

  if (npoint == 0) error->all(FLERR,"Surf file does not contain points");
  if (dim == 2 && nline == 0) 
    error->all(FLERR,"Surf file does not contain lines");
  if (dim == 3 && ntri == 0) 
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

  int nread = 0;
  
  while (nread < npoint) {
    if (npoint-nread > CHUNK) nchunk = CHUNK;
    else nchunk = npoint-nread;
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
      pts[npoint].x[0] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
      pts[npoint].x[1] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
      if (dim == 3) 
        pts[npoint].x[2] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
      else pts[npoint].x[2] = 0.0;
      npoint++;
      buf = next + 1;
    }

    nread += nchunk;
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  %d points\n",npoint);
    if (logfile) fprintf(logfile,"  %d points\n",npoint);
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

  int nread = 0;

  while (nread < nline) {
    if (nline-nread > CHUNK) nchunk = CHUNK;
    else nchunk = nline-nread;
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

    // allow for optional type in each line element

    if (nwords != 3 && nwords != 4)
      error->all(FLERR,"Incorrect line format in surf file");
    int typeflag = 0;
    if (nwords == 4) typeflag = 1;

    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');
      strtok(buf," \t\n\r\f");
      if (typeflag) type = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
      else type = 1;
      p1 = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
      p2 = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
      if (p1 < 1 || p1 > npoint || p2 < 1 || p2 > npoint || p1 == p2)
	error->all(FLERR,"Invalid point index in line");
      lines[nline].type = type;
      lines[nline].mask = 1;
      lines[nline].p1 = p1-1;
      lines[nline].p2 = p2-1;
      nline++;
      buf = next + 1;
    }

    nread += nchunk;
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  %d lines\n",nline);
    if (logfile) fprintf(logfile,"  %d lines\n",nline);
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

  int nread = 0;

  while (nread < ntri) {
    if (ntri-nread > CHUNK) nchunk = CHUNK;
    else nchunk = ntri-nread;
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

    // allow for optional type in each tri element

    if (nwords != 4 && nwords != 5)
      error->all(FLERR,"Incorrect line format in surf file");
    int typeflag = 0;
    if (nwords == 5) typeflag = 1;

    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');
      strtok(buf," \t\n\r\f");
      if (typeflag) type = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
      else type = 1;
      p1 = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
      p2 = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
      p3 = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
      if (p1 < 1 || p1 > npoint || p2 < 1 || p2 > npoint || 
	  p3 < 1 || p3 > npoint || p1 == p2 || p2 == p3)
	error->all(FLERR,"Invalid point index in triangle");
      tris[ntri].type = type;
      tris[ntri].mask = 1;
      tris[ntri].p1 = p1-1;
      tris[ntri].p2 = p2-1;
      tris[ntri].p3 = p3-1;
      ntri++;
      buf = next + 1;
    }

    nread += nchunk;
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  %d triangles\n",ntri);
    if (logfile) fprintf(logfile,"  %d triangles\n",ntri);
  }
}

/* ----------------------------------------------------------------------
   translate new vertices by (dx,dy,dz)
   for 2d, dz will be 0.0
------------------------------------------------------------------------- */

void ReadSurf::translate(double dx, double dy, double dz)
{
  for (int i = 0; i < npoint; i++) {
    pts[i].x[0] += dx;
    pts[i].x[1] += dy;
    pts[i].x[2] += dz;
  }
}

/* ----------------------------------------------------------------------
   scale new vertices by (sx,sy,sz) around origin
   for 2d, do not reset x[2] to avoid epsilon change
------------------------------------------------------------------------- */

void ReadSurf::scale(double sx, double sy, double sz)
{
  for (int i = 0; i < npoint; i++) {
    pts[i].x[0] = sx*(pts[i].x[0]-origin[0]) + origin[0];
    pts[i].x[1] = sy*(pts[i].x[1]-origin[1]) + origin[1];
    if (dim == 3) pts[i].x[2] = sz*(pts[i].x[2]-origin[2]) + origin[2];
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

  for (int i = 0; i < npoint; i++) {
    d[0] = pts[i].x[0] - origin[0];
    d[1] = pts[i].x[1] - origin[1];
    d[2] = pts[i].x[2] - origin[2];
    MathExtra::matvec(rotmat,d,dnew);
    pts[i].x[0] = dnew[0] + origin[0];
    pts[i].x[1] = dnew[1] + origin[1];
    if (dim == 3) pts[i].x[2] = dnew[2] + origin[2];
  }
}

/* ----------------------------------------------------------------------
   invert new vertex ordering within each line or tri
   this flips direction of surface normal
------------------------------------------------------------------------- */

void ReadSurf::invert()
{
  int tmp;

  if (dim == 2) {
    int m = nline_old;
    for (int i = 0; i < nline_new; i++) {
      tmp = lines[m].p1;
      lines[m].p1 = lines[m].p2;
      lines[m].p2 = tmp;
      m++;
    }
  }

  if (dim == 3) {
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
  int i,n,nptotal,nltotal,dim,side,flag;
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
    // nptotal,nltotal = total current count of points and lines, including old

    nptotal = npoint;
    nltotal = nline;

    for (i = 0; i < nline; i++) {
      x1 = pts[lines[i].p1].x;
      x2 = pts[lines[i].p2].x;
      flag = 0;
      if (x1[dim] < value && x2[dim] > value) flag = 1;
      if (x1[dim] > value && x2[dim] < value) flag = 1;
      if (flag) {
	if (nptotal == maxpoint) {
	  maxpoint += DELTA;
	  pts = (Point *) 
	    memory->srealloc(pts,maxpoint*sizeof(Point),"readsurf:pts");
	  x1 = pts[lines[i].p1].x;
	  x2 = pts[lines[i].p2].x;
	}
	if (nltotal == maxline) {
	  maxline += DELTA;
	  lines = (Line *) 
	    memory->srealloc(lines,maxline*sizeof(Line),"readsurf:lines");
	}

	param = (value-x1[dim]) / (x2[dim]-x1[dim]);
	pts[nptotal].x[0] = x1[0] + param*(x2[0]-x1[0]);
	pts[nptotal].x[1] = x1[1] + param*(x2[1]-x1[1]);
	pts[nptotal].x[2] = 0.0;
	pts[nptotal].x[dim] = value;

	memcpy(&lines[nltotal],&lines[i],sizeof(Surf::Line));
	lines[i].p2 = nptotal;
	lines[nltotal].p1 = nptotal;

        nptotal++;
        nltotal++;
      }
    }

    npoint = nptotal;
    nline = nltotal;

    // project all points outside clipping edge to the edge

    for (i = 0; i < npoint; i++) {
      x = pts[i].x;
      if (side == 0 && x[dim] < value) x[dim] = value;
      if (side == 1 && x[dim] > value) x[dim] = value;
    }

    // ptflag[I] = # of lines that include point I
    
    int *ptflag;
    memory->create(ptflag,npoint,"readsurf:ptflag");
    for (i = 0; i < npoint; i++) ptflag[i] = 0;

    for (i = 0; i < nline; i++) {
      ptflag[lines[i].p1]++;
      ptflag[lines[i].p2]++;
    }

    // lineflag[I] = 1 if line I can be removed b/c lies on clip edge
    // also decrement ptflag for pts in removed line
    
    int *lineflag;
    memory->create(lineflag,nline,"readsurf:lineflag");
    for (i = 0; i < nline; i++) lineflag[i] = 0;

    for (i = 0; i < nline; i++) {
      if (pts[lines[i].p1].x[dim] == value && 
          pts[lines[i].p2].x[dim] == value) {
	lineflag[i] = 1;
	ptflag[lines[i].p1]--;
	ptflag[lines[i].p2]--;
      }
    }

    // condense pts/lines to remove deleted points and lines
    // when delete points, set ptflag[I] = new index of point I for kept points
    // renumber point indices in lines using altered ptflag
    // reset npoint_new,nline_new

    n = 0;
    for (i = 0; i < npoint; i++) {
      if (ptflag[i]) {
	memcpy(&pts[n],&pts[i],sizeof(Point));
	ptflag[i] = n;
	n++;
      }
    }
    npoint = n;
    
    n = 0;
    for (i = 0; i < nline; i++) {
      if (!lineflag[i]) {
	memcpy(&lines[n],&lines[i],sizeof(Line));
	lines[n].p1 = ptflag[lines[n].p1];
	lines[n].p2 = ptflag[lines[n].p2];
	n++;
      }
    }
    nline = n;

    // clean up

    memory->destroy(ptflag);
    memory->destroy(lineflag);
  }

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  clipped to %d points\n",npoint);
      fprintf(screen,"  clipped to %d lines\n",nline);
    }
    if (logfile) {
      fprintf(logfile,"  clipped to %d points\n",npoint);
      fprintf(logfile,"  clipped to %d lines\n",nline);
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
  int i,n,nptotal,nttotal,dim,side;
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
    // nptotal,nttotal = total current count of points and tris, including old

    nedge = maxedge = 0;
    edge = NULL;
    
    nptotal = npoint;
    nttotal = ntri;

    for (i = 0; i < ntri; i++) {
      p[0] = tris[i].p1;
      p[1] = tris[i].p2;
      p[2] = tris[i].p3;
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
	if (nptotal+2 >= maxpoint) {
	  maxpoint += DELTA;
	  pts = (Point *) 
	    memory->srealloc(pts,maxpoint*sizeof(Point),"surf:pts");
	  xp[0] = pts[p[0]].x;
	  xp[1] = pts[p[1]].x;
	  xp[2] = pts[p[2]].x;
	}
	if (nttotal == maxtri) {
	  maxtri += DELTA;
	  tris = (Tri *) 
	    memory->srealloc(tris,maxtri*sizeof(Tri),"surf:tris");
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
	    pts[nptotal].x[0] = xp[good][0] + param*(xp[bad][0]-xp[good][0]);
	    pts[nptotal].x[1] = xp[good][1] + param*(xp[bad][1]-xp[good][1]);
	    pts[nptotal].x[2] = xp[good][2] + param*(xp[bad][2]-xp[good][2]);
	    pts[nptotal].x[dim] = value;
	    ipt = nptotal;
	    add_edge(p[good],p[bad],ipt);
            nptotal++;
	  }

	  if (bad == 0) tris[i].p1 = ipt;
	  else if (bad == 1) tris[i].p2 = ipt;
	  else if (bad == 2) tris[i].p3 = ipt;

	  if (nbad == 2) {
	    if (bad == 0 && pflag[1] == BAD) bad = 1;
	    else bad = 2;

	    ipt = find_edge(p[good],p[bad]);
	    if (ipt < 0) {
	      param = (value-xp[good][dim]) / (xp[bad][dim]-xp[good][dim]);
	      pts[nptotal].x[0] = 
		xp[good][0] + param*(xp[bad][0]-xp[good][0]);
	      pts[nptotal].x[1] = 
		xp[good][1] + param*(xp[bad][1]-xp[good][1]);
	      pts[nptotal].x[2] = 
		xp[good][2] + param*(xp[bad][2]-xp[good][2]);
	      pts[nptotal].x[dim] = value;
	      ipt = nptotal;
	      add_edge(p[good],p[bad],nptotal);
              nptotal++;
	    }

	    if (bad == 0) tris[i].p1 = ipt;
	    else if (bad == 1) tris[i].p2 = ipt;
	    else if (bad == 2) tris[i].p3 = ipt;
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
	    pts[nptotal].x[0] = xp[good][0] + param*(xp[bad][0]-xp[good][0]);
	    pts[nptotal].x[1] = xp[good][1] + param*(xp[bad][1]-xp[good][1]);
	    pts[nptotal].x[2] = xp[good][2] + param*(xp[bad][2]-xp[good][2]);
	    pts[nptotal].x[dim] = value;
	    ipt = nptotal;
	    add_edge(p[good],p[bad],ipt);
            nptotal++;
	  }

	  if (bad == 0) tris[i].p1 = ipt;
	  else if (bad == 1) tris[i].p2 = ipt;
	  else if (bad == 2) tris[i].p3 = ipt;

	  memcpy(&tris[nttotal],&tris[i],sizeof(Tri));
	  if (good == 0) tris[nttotal].p1 = ipt;
	  else if (good == 1) tris[nttotal].p2 = ipt;
	  else if (good == 2) tris[nttotal].p3 = ipt;

	  if (good == 0 && pflag[1] == GOOD) good = 1;
	  else good = 2;

	  ipt = find_edge(p[good],p[bad]);
	  if (ipt < 0) {
	    param = (value-xp[good][dim]) / (xp[bad][dim]-xp[good][dim]);
	    pts[nptotal].x[0] = xp[good][0] + param*(xp[bad][0]-xp[good][0]);
	    pts[nptotal].x[1] = xp[good][1] + param*(xp[bad][1]-xp[good][1]);
	    pts[nptotal].x[2] = xp[good][2] + param*(xp[bad][2]-xp[good][2]);
	    pts[nptotal].x[dim] = value;
	    ipt = nptotal;
	    add_edge(p[good],p[bad],ipt);
            nptotal++;
	  }

	  if (bad == 0) tris[nttotal].p1 = ipt;
	  else if (bad == 1) tris[nttotal].p2 = ipt;
	  else if (bad == 2) tris[nttotal].p3 = ipt;
	  nttotal++;
	}
      }
    }

    memory->destroy(edge);

    npoint = nptotal;
    ntri = nttotal;

    // project all points outside clipping plane to the plane

    for (i = 0; i < npoint; i++) {
      x = pts[i].x;
      if (side == 0 && x[dim] < value) x[dim] = value;
      if (side == 1 && x[dim] > value) x[dim] = value;
    }

    // ptflag[I] = # of tris that include point I
    
    int *ptflag;
    memory->create(ptflag,npoint,"readsurf:ptflag");
    for (i = 0; i < npoint; i++) ptflag[i] = 0;

    for (i = 0; i < ntri; i++) {
      ptflag[tris[i].p1]++;
      ptflag[tris[i].p2]++;
      ptflag[tris[i].p3]++;
    }

    // triflag[I] = 1 if tri I can be removed b/c lies on clip plane
    // also decrement ptflag for pts in removed tri
    
    int *triflag;
    memory->create(triflag,ntri,"readsurf:triflag");
    for (i = 0; i < ntri; i++) triflag[i] = 0;

    for (i = 0; i < ntri; i++) {
      if (pts[tris[i].p1].x[dim] == value && 
	  pts[tris[i].p2].x[dim] == value &&
	  pts[tris[i].p3].x[dim] == value) {
	triflag[i] = 1;
	ptflag[tris[i].p1]--;
	ptflag[tris[i].p2]--;
	ptflag[tris[i].p3]--;
      }
    }

    // condense pts/tris to remove deleted points and tris
    // when delete points, set ptflag[I] = new index of point I for kept points
    // renumber point indices in tris using altered ptflag
    // reset npoint_new,ntri_new

    n = 0;
    for (i = 0; i < npoint; i++) {
      if (ptflag[i]) {
	memcpy(&pts[n],&pts[i],sizeof(Point));
	ptflag[i] = n;
	n++;
      }
    }
    npoint = n;
    
    n = 0;
    for (i = 0; i < ntri; i++) {
      if (!triflag[0]) {
	memcpy(&tris[n],&tris[i],sizeof(Tri));
	tris[n].p1 = ptflag[tris[n].p1];
	tris[n].p2 = ptflag[tris[n].p2];
	tris[n].p3 = ptflag[tris[n].p3];
	n++;
      }
    }
    ntri = n;

    // clean up

    memory->destroy(ptflag);
    memory->destroy(triflag);
  }

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  clipped to %d points\n",npoint);
      fprintf(screen,"  clipped to %d tris\n",ntri);
    }
    if (logfile) {
      fprintf(logfile,"  clipped to %d points\n",npoint);
      fprintf(logfile,"  clipped to %d tris\n",ntri);
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

  for (int i = 0; i < npoint; i++) {
    x = pts[i].x;
    if (x[0] >= boxlo[0] && x[0] <= boxhi[0]) {
      if (x[0]-boxlo[0] < xdelta) x[0] = boxlo[0];
      else if (boxhi[0]-x[0] < xdelta) x[0] = boxhi[0];
    }
    if (x[1] >= boxlo[1] && x[1] <= boxhi[1]) {
      if (x[1]-boxlo[1] < ydelta) x[1] = boxlo[1];
      else if (boxhi[1]-x[1] < ydelta) x[1] = boxhi[1];
    }
    if (dim == 2) continue;
    if (x[2] >= boxlo[2] && x[2] <= boxhi[2]) {
      if (x[2]-boxlo[2] < zdelta) x[2] = boxlo[2];
      else if (boxhi[2]-x[2] < zdelta) x[2] = boxhi[2];
    }
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
  memory->create(count,npoint,"readsurf:count");
  memory->create(p2e,npoint,2,"readsurf:count");
  for (int i = 0; i < npoint; i++) count[i] = 0;

  for (int i = 0; i < nline; i++) {
    p1 = lines[i].p1;
    p2 = lines[i].p2;
    p2e[p1][count[p1]++] = i;
    p2e[p2][count[p2]++] = i;
  }
  
  // check that norms of adjacent lines are not in opposite directions
  // norms are stored in Surf::lines, at end of orignal nline_old surfs

  Surf::Line *surflines = surf->lines;

  double dot;
  double *norm1,*norm2;

  int nerror = 0;
  int nwarn = 0;
  for (int i = 0; i < npoint; i++) {
    if (count[i] == 1) continue;
    norm1 = surflines[p2e[i][0] + nline_old].norm;
    norm2 = surflines[p2e[i][1] + nline_old].norm;
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
#elif defined SPARTA_UNORDERED_MAP
  std::unordered_map<bigint,int> hash;
  std::unordered_map<bigint,int>::iterator it;
#else
  std::tr1::unordered_map<bigint,int> hash;
  std::tr1::unordered_map<bigint,int>::iterator it;
#endif

  // insert each edge into hash with triangle as value

  bigint p1,p2,p3,key;

  for (int i = 0; i < ntri; i++) {
    p1 = tris[i].p1;
    p2 = tris[i].p2;
    p3 = tris[i].p3;
    key = (p1 << 32) | p2;
    hash[key] = i;
    key = (p2 << 32) | p3;
    hash[key] = i;
    key = (p3 << 32) | p1;
    hash[key] = i;
  }

  // check that norms of adjacent triangles are not in opposite directions
  // norms are stored in Surf::tris, at end of orignal ntri_old surfs

  Surf::Tri *surftris = surf->tris;

  double dot;
  double *norm1,*norm2;

  int nerror = 0;
  int nwarn = 0;
  for (it = hash.begin(); it != hash.end(); ++it) {
    p1 = it->first >> 32;
    p2 = it->first & MAXSMALLINT;
    key = (p2 << 32) | p1;
    if (hash.find(key) == hash.end()) continue;
    norm1 = surftris[it->second + ntri_old].norm;
    norm2 = surftris[hash[key] + ntri_old].norm;
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
  int i,j,n;
  int *csurfs;
  double side,epssq;
  double *p1,*p2,*lo,*hi;
  Surf::Line *line;

  Surf::Line *surflines = surf->lines;
  Grid::ChildCell *cells = grid->cells;
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
      line = &surflines[csurfs[i]];
      for (j = 0; j < n; j++) {
        if (i == j) continue;
        p1 = surflines[csurfs[j]].p1;
        p2 = surflines[csurfs[j]].p2;
        point_line_compare(p1,line->p1,line->p2,epssq,nerror,nwarn);
        point_line_compare(p2,line->p1,line->p2,epssq,nerror,nwarn);
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
  int i,j,n;
  int *csurfs;
  double side,epssq;
  double *p1,*p2,*p3,*lo,*hi;
  Surf::Tri *tri;

  Surf::Tri *surftris = surf->tris;
  Grid::ChildCell *cells = grid->cells;
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
      tri = &surftris[csurfs[i]];
      for (j = 0; j < n; j++) {
        if (i == j) continue;
        p1 = surftris[csurfs[j]].p1;
        p2 = surftris[csurfs[j]].p2;
        p3 = surftris[csurfs[j]].p3;
        point_tri_compare(p1,tri->p1,tri->p2,tri->p3,tri->norm,
                          epssq,nerror,nwarn,icell,csurfs[i],csurfs[j]);
        point_tri_compare(p2,tri->p1,tri->p2,tri->p3,tri->norm,
                          epssq,nerror,nwarn,icell,csurfs[i],csurfs[j]);
        point_tri_compare(p3,tri->p1,tri->p2,tri->p3,tri->norm,
                          epssq,nerror,nwarn,icell,csurfs[i],csurfs[j]);
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

void ReadSurf::point_line_compare(double *pt, double *p1, double *p2,
                                  double epssq, int &nerror, int &nwarn)
{
  if (pt[0] == p1[0] && pt[1] == p1[1]) return;
  if (pt[0] == p2[0] && pt[1] == p2[1]) return;
  double rsq = Geometry::distsq_point_line(pt,p1,p2);
  if (rsq == 0.0) nerror++;
  else if (rsq < epssq) nwarn++;
}

/* ----------------------------------------------------------------------
   compute distance bewteen a point and triangle
   just return if point is an endpoint of triangle
   increment nerror if point on triangle
   increment nwarn if point is within epssq distance of triangle
------------------------------------------------------------------------- */

void ReadSurf::point_tri_compare(double *pt, double *p1, double *p2, double *p3,
                                 double *norm,
                                 double epssq, int &nerror, int &nwarn,
                                 int, int, int) 
{
  if (pt[0] == p1[0] && pt[1] == p1[1] && pt[2] == p1[2]) return;
  if (pt[0] == p2[0] && pt[1] == p2[1] && pt[2] == p2[2]) return;
  if (pt[0] == p3[0] && pt[1] == p3[1] && pt[2] == p3[2]) return;
  double rsq = Geometry::distsq_point_tri(pt,p1,p2,p3,norm);
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
  if (dim == 3) epsilon = MIN(epsilon,domain->zprd);
  epsilon *= EPSILON;
  double epssq = epsilon * epsilon;

  // goal: N roughly cubic bins where N = # of new points
  // nbinxyz = # of bins in each dim
  // xyzbin = bin size in each dim
  // for 2d, nbinz = 1
  // after setting bin size, add 1 to nbinxyz
  // this allows for 2nd binning via offset origin

  int nbinx,nbiny,nbinz;
  double xbin,ybin,zbin;
  double xbininv,ybininv,zbininv;
  
  if (dim == 2) {
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
  // do not offset bin origin in a dim with only 1 bin

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
