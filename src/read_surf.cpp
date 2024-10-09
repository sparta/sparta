/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
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
#include "dirent.h"
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
#define CHUNK 16384
#define EPSILON_NORM 1.0e-12
#define BIG 1.0e20
#define DELTA 1024
#define DELTA_DISCARD 128

/* ---------------------------------------------------------------------- */

ReadSurf::ReadSurf(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

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
  if (surf->implicit)
    error->all(FLERR,"Cannot use read_surf when global surfs implicit is set");

  surf->exist = 1;
  dim = domain->dimension;
  distributed = surf->distributed;

  if (narg < 1) error->all(FLERR,"Illegal read_surf command");

  // if filename contains "*", search dir for latest surface file

  char *file = new char[strlen(arg[0]) + 16];
  if (strchr(arg[0],'*')) {
    int n;
    if (me == 0) {
      file_search(arg[0],file);
      n = strlen(file) + 1;
    }
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(file,n,MPI_CHAR,0,world);
  } else strcpy(file,arg[0]);

  // check for multiproc files

  if (strchr(arg[0],'%')) multiproc = 1;
  else multiproc = 0;

  // check for initial "type" and "custom" keywords
  // they affect file format

  int iarg = 1;

  typeflag = 0;
  ncustom = 0;
  name_custom = NULL;
  type_custom = NULL;
  size_custom = NULL;
  index_custom = NULL;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"type") == 0) {
      typeflag = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"custom") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Invalid read_surf command");

      name_custom = (char **)
	memory->srealloc(name_custom,(ncustom+1)*sizeof(char *),
			 "readsurf:name_custom");
      memory->grow(type_custom,ncustom+1,"readsurf:type_custom");
      memory->grow(size_custom,ncustom+1,"readsurf:size_custom");
      memory->grow(index_custom,ncustom+1,"readsurf:index_custom");

      int n = strlen(arg[iarg+1]) + 1;
      name_custom[ncustom] = new char[n];
      strcpy(name_custom[ncustom],arg[iarg+1]);
      if (strcmp(arg[iarg+2],"int") == 0) type_custom[ncustom] = 0;
      else if (strcmp(arg[iarg+2],"float") == 0) type_custom[ncustom] = 1;
      else error->all(FLERR,"Invalid read_surf command");
      size_custom[ncustom] = input->inumeric(FLERR,arg[iarg+3]);
      if (size_custom[ncustom] < 0)
	error->all(FLERR,"Invalid read_surf command");
      ncustom++;

      iarg += 4;
    } else break;
  }

  // setup for temporary storage of per-surf custom values in file
  // nvalues_custom = # of custom values per surf
  // cvalues = growable array with ID + custom values for each surf
  // if custom name already exists, check compatibility of type/size args

  nvalues_custom = 0;
  cvalues = NULL;

  if (ncustom) {
    for (int ic = 0; ic < ncustom; ic++) {
      int index = surf->find_custom(name_custom[ic]);
      if (index >= 0) {
	int flag = 0;
	if (type_custom[ic] != surf->etype[index]) flag = 1;
	if (size_custom[ic] != surf->esize[index]) flag = 1;
	if (flag) error->all(FLERR,"Read_surf custom attributes do not match "
	 		     "already existing custom data");
	index_custom[ic] = index;
        surf->estatus[index] = 0;
      } else {
	index_custom[ic] =
	  surf->add_custom(name_custom[ic],type_custom[ic],size_custom[ic]);
      }

      if (size_custom[ic] == 0) nvalues_custom++;
      else nvalues_custom += size_custom[ic];
    }
  }

  // -----------------------
  // read surface data from file(s)
  // -----------------------

  if (me == 0)
    if (screen) fprintf(screen,"Reading surface file ...\n");

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  // multiproc = 0/1 = single or multiple files
  // each file may list Points or not
  // store surfs in distributed fashion in local lines/tris for now

  nsurf = maxsurf = 0;
  lines = NULL;
  tris = NULL;

  if (!multiproc) read_single(file);
  else read_multiple(file);

  delete [] file;

  // print stats

  if (me == 0) {
    if (dim == 2) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " lines\n",nsurf_all);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " lines\n",nsurf_all);
    }
    if (dim == 3) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " triangles\n",nsurf_all);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " triangles\n",nsurf_all);
    }
  }

  // check for surf IDs within bounds

  check_bounds();

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // -----------------------
  // transform and check surface elements
  // -----------------------

  // parse command-line args after type/custom
  // apply geometric transformations when parsed
  // apply keywords group,typeadd,transparent,particle,file below

  process_args(iarg,narg,arg);

  // add surface elements to specified group

  if (grouparg) {
    int igroup = surf->find_group(arg[grouparg]);
    if (igroup < 0) igroup = surf->add_group(arg[grouparg]);
    int groupbit = surf->bitmask[igroup];
    if (dim == 2)
      for (int i = 0; i < nsurf; i++) lines[i].mask |= groupbit;
    else
      for (int i = 0; i < nsurf; i++) tris[i].mask |= groupbit;
  }

  // offset surface types

  if (typeadd) {
    if (dim == 2)
      for (int i = 0; i < nsurf; i++) lines[i].type += typeadd;
    else
      for (int i = 0; i < nsurf; i++) tris[i].type += typeadd;
  }

  // flag surfs as transparent

  if (transparent_flag) {
    if (dim == 2)
      for (int i = 0; i < nsurf; i++) lines[i].transparent = 1;
    else
      for (int i = 0; i < nsurf; i++) tris[i].transparent = 1;
  }

  // list of read-in surfs in lines/tris/cvalues is now final (post-clip)
  // add all read-in data to Surf data structs via add_surfs()
  // can then deallocate all local data

  bigint nsurf_old = surf->nsurf;
  int nsurf_old_mine;
  if (!distributed) nsurf_old_mine = surf->nlocal;
  else nsurf_old_mine = surf->nown;

  surf->add_surfs(0,nsurf,lines,tris,ncustom,index_custom,cvalues);

  memory->sfree(lines);
  memory->sfree(tris);

  if (ncustom) {
    for (int ic = 0; ic < ncustom; ic++) delete [] name_custom[ic];
    memory->sfree(name_custom);
    memory->destroy(type_custom);
    memory->destroy(size_custom);
    memory->destroy(index_custom);
    memory->destroy(cvalues);
  }

  // write out new surf file if requested
  // do this before grid cell assignment or checks, in case error occurs

  if (filearg) {
    WriteSurf *wf = new WriteSurf(sparta);
    wf->statflag = 0;
    wf->command(narg-filearg,&arg[filearg]);
    delete wf;
  }

  // output extent of new surfs, tiny ones may have been created by clip

  if (dim == 2) surf->output_extent(nsurf_old_mine);
  else surf->output_extent(nsurf_old_mine);

  // compute normals of new surfs

  if (dim == 2) surf->compute_line_normal(nsurf_old_mine);
  else surf->compute_tri_normal(nsurf_old_mine);

  // error check on new surfs
  // all points must be inside or on surface of simulation box

  if (dim == 2) surf->check_point_inside(nsurf_old_mine);
  else surf->check_point_inside(nsurf_old_mine);

  // error checks that can be done before surfs are mapped to grid cells
  // need to do watertight test on old + new surfs
  // could do neighbor norm test only on new surfs (currently a no-op)

  if (dim == 2) {
    surf->check_watertight_2d();
    check_neighbor_norm_2d();
  } else {
    surf->check_watertight_3d();
    check_neighbor_norm_3d();
  }

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // -----------------------
  // map surfs to grid cells
  // -----------------------

  // sort particles

  if (particle->exist) particle->sort();

  // make list of surf elements I own
  // clear grid of surf info including split cells

  surf->setup_owned();
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
  double time4 = MPI_Wtime();

  // assign surfs to grid cells

  grid->surf2grid(1);

  MPI_Barrier(world);
  double time5 = MPI_Wtime();

  // error check on any points too near other surfs
  // done on per-grid-cell basis, expensive to do globally

  if (dim == 2) surf->check_point_near_surf_2d();
  else surf->check_point_near_surf_3d();

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
        surfint *csurfs = cells[icell].csurfs;
        int m;
        if (dim == 2) {
          Surf::Line *lines = surf->lines;
          for (m = 0; m < nsurf; m++) {
            if (lines[csurfs[m]].id > nsurf_old) break;
          }
        } else {
          Surf::Tri *tris = surf->tris;
          for (m = 0; m < nsurf; m++) {
            if (tris[csurfs[m]].id >= nsurf_old) break;
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

  // stats

  double time_total = time6-time1;
  double time_s2g = time5-time4;

  if (comm->me == 0) {
    if (screen) {
      if (particle->exist)
        fprintf(screen,"  " BIGINT_FORMAT " deleted particles\n",ndeleted);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  read/check/sort/surf2grid/ghost/"
              "inout/particle percent = "
              "%g %g %g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total,100.0*(time7-time6)/time_total,
              100.0*(time8-time7)/time_total);
      fprintf(screen,"  surf2grid time = %g secs\n",time_s2g);
      fprintf(screen,"  map/comm1/comm2/comm3/comm4/split percent = "
              "%g %g %g %g %g %g\n",
              100.0*grid->tmap/time_s2g,100.0*grid->tcomm1/time_s2g,
              100.0*grid->tcomm2/time_s2g,100.0*grid->tcomm3/time_s2g,
              100.0*grid->tcomm4/time_s2g,100.0*grid->tsplit/time_s2g);
    }

    if (logfile) {
      if (particle->exist)
        fprintf(logfile,"  " BIGINT_FORMAT " deleted particles\n",ndeleted);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  read/check/sort/surf2grid/ghost/"
              "inout/particle percent = "
              "%g %g %g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total,100.0*(time7-time6)/time_total,
              100.0*(time8-time7)/time_total);
      fprintf(logfile,"  surf2grid time = %g secs\n",time_s2g);
      fprintf(logfile,"  map/comm1/comm2/comm3/comm4/split percent = "
              "%g %g %g %g %g %g\n",
              100.0*grid->tmap/time_s2g,100.0*grid->tcomm1/time_s2g,
              100.0*grid->tcomm2/time_s2g,100.0*grid->tcomm3/time_s2g,
              100.0*grid->tcomm4/time_s2g,100.0*grid->tsplit/time_s2g);
    }
  }
}

// ----------------------------------------------------------------------
// methods to read surface files
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   read a single input file
   store surfs directly in lines/tris (all) or mylines/mytris (distributed)
------------------------------------------------------------------------- */

void ReadSurf::read_single(char *file)
{
  // read file

  if (me == 0) filereader = 1;
  else filereader = 0;
  filecomm = world;
  me_file = me;
  nprocs_file = nprocs;

  read_file(file);

  nsurf_all = nsurf_file;
}

/* ----------------------------------------------------------------------
   read multiple files, cluster of procs for each file
   store surfs initially in distributed tmplines/tmptris
   communicate to populate lines/tris (all) or mylines/mytris (distributed)
------------------------------------------------------------------------- */

void ReadSurf::read_multiple(char *file)
{
  // read base file, sets nsurf_all

  base(file);

  // make procfile copy so can substitute for "%"

  char *procfile = new char[strlen(file) + 16];
  char *ptr = strchr(file,'%');

  // if nprocs > files, break procs into clusters, each cluster reads one file

  bigint nsurf_allfiles;

  if (nprocs > nfiles) {

    // icluster = which cluster of procs this proc is in
    // fileproc = ID of proc in my cluster who reads from file
    // filereader = 1 if this proc reads file, else 0

    int icluster = static_cast<int> ((bigint) me * nfiles/nprocs);
    int fileproc = static_cast<int> ((bigint) icluster * nprocs/nfiles);
    int fcluster = static_cast<int> ((bigint) fileproc * nfiles/nprocs);
    if (fcluster < icluster) fileproc++;
    if (me == fileproc) filereader = 1;
    else filereader = 0;
    MPI_Comm_split(world,icluster,0,&filecomm);
    MPI_Comm_rank(filecomm,&me_file);
    MPI_Comm_size(filecomm,&nprocs_file);

    // each cluster reads a file, stores surfs in tmplines/tmptris

    *ptr = '\0';
    sprintf(procfile,"%s%d%s",file,icluster,ptr+1);
    *ptr = '%';

    read_file(procfile);

    bigint nsurf_onefile = 0;
    if (filereader) nsurf_onefile = nsurf_file;
    MPI_Allreduce(&nsurf_onefile,&nsurf_allfiles,1,MPI_SPARTA_BIGINT,
		  MPI_SUM,world);

  // if nprocs <= files, each proc reads one or more files

  } else {

    filereader = 1;
    MPI_Comm_split(world,me,0,&filecomm);
    me_file = 0;
    nprocs_file = 1;

    // each proc reads every Pth file, stores surfs in tmplines/tmptris

    bigint nsurf_oneproc = 0;

    for (int iproc = me; iproc < nfiles; iproc += nprocs) {
      *ptr = '\0';
      sprintf(procfile,"%s%d%s",file,iproc,ptr+1);
      *ptr = '%';

      read_file(procfile);
      nsurf_oneproc += nsurf_file;
    }

    MPI_Allreduce(&nsurf_oneproc,&nsurf_allfiles,1,MPI_SPARTA_BIGINT,
		  MPI_SUM,world);
  }

  delete [] procfile;

  // check that read-in surf count matches base file

  if (nsurf_allfiles != nsurf_all) {
    char str[128];
    fprintf(screen,"Read surf mismatch in surf count across mutiple files: base "
	    BIGINT_FORMAT ", actual " BIGINT_FORMAT "\n",nsurf_all,nsurf_allfiles);
    error->all(FLERR,str);
  }

  // clean-up

  MPI_Comm_free(&filecomm);
}

/* ----------------------------------------------------------------------
   read a single input file, header and post-header sections
   proc with filereader = 1 reads the file
   filecomm = MPI communicator for all procs who share this file
   store surfs in lines/tris
------------------------------------------------------------------------- */

void ReadSurf::read_file(char *file)
{
  // open file

  if (filereader) open(file);

  // read header

  pts = NULL;
  header();

  // read and store data from Points section

  parse_keyword(1);

  if (strcmp(keyword,"Points") == 0) {
    if (npoint_file == 0)
      error->all(FLERR,"Read_surf file has no points keyword");
    read_points();
    parse_keyword(0);
  } else if (npoint_file)
    error->all(FLERR,"Read_surf file has no Points section");

  // read and store data from Lines or Triangles section

  if (dim == 2) {
    if (strcmp(keyword,"Lines") != 0)
      error->all(FLERR,"Read_surf did not find Lines section of surf file");
    read_lines();
  } else {
    if (strcmp(keyword,"Triangles") != 0)
      error->all(FLERR,"Read_surf did not find Triangles section of surf file");
    read_tris();
  }

  // can now free Points

  if (npoint_file) memory->sfree(pts);

  // close file

  if (filereader) {
    if (compressed) pclose(fp);
    else fclose(fp);
  }
}

/* ----------------------------------------------------------------------
   read free-format base file for multiproc files
   1st line and blank lines are skipped
   non-blank lines are checked for header keywords and leading value is read
------------------------------------------------------------------------- */

void ReadSurf::base(char *file)
{
  int n;
  char *ptr;

  // open base file

  if (me == 0) {
    char *hfile;
    hfile = new char[strlen(file) + 16];
    char *ptr = strchr(file,'%');
    *ptr = '\0';
    sprintf(hfile,"%s%s%s",file,"base",ptr+1);
    *ptr = '%';
    fp = fopen(hfile,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open surface base file %s",hfile);
      error->one(FLERR,str);
    }
    delete [] hfile;
  }

  // skip 1st line of file

  if (me == 0) {
    char *eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of surf base file");
  }

  // read keywords from file

  nfiles = 0;
  bigint nline_basefile = 0;
  bigint ntri_basefile = 0;

  while (1) {

    // read a line and bcast length

    if (me == 0) {
      if (fgets(line,MAXLINE,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
    }
    MPI_Bcast(&n,1,MPI_INT,0,world);

    // if n = 0 then end-of-file so break

    if (n == 0) break;

    // bcast line

    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // trim anything from '#' onward
    // if line is blank, continue

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    // search line for header keyword and set corresponding variable

    if (strstr(line,"files")) {
      sscanf(line,"%d",&nfiles);
    } else if (strstr(line,"lines")) {
      if (dim == 3)
        error->all(FLERR,"Surf file cannot contain lines for 3d simulation");
      sscanf(line,BIGINT_FORMAT,&nline_basefile);
      nsurf_all = nline_basefile;
    } else if (strstr(line,"triangles")) {
      if (dim == 2)
        error->all(FLERR,
                   "Surf file cannot contain triangles for 2d simulation");
      sscanf(line,BIGINT_FORMAT,&ntri_basefile);
      nsurf_all = ntri_basefile;

    } else error->all(FLERR,"Invalid keyword in surf base file");
  }

  if (nfiles == 0)
    error->all(FLERR,"Surf base file does not contain files");
  if (dim == 2 && nline_basefile == 0)
    error->all(FLERR,"Surf base file does not contain lines");
  if (dim == 3 && ntri_basefile == 0)
    error->all(FLERR,"Surf base file does not contain triangles");

  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   read free-format header of data file
   1st line and blank lines are skipped
   non-blank lines are checked for header keywords and leading value is read
   header ends with non-blank line containing no header keyword (or EOF)
   return line with non-blank line (or empty line if EOF)
   file does not have to contain points keyword
------------------------------------------------------------------------- */

void ReadSurf::header()
{
  int n;
  char *ptr;

  // skip 1st line of file

  if (filereader) {
    char *eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of surf file");
  }

  npoint_file = 0;
  int nline_file = 0;
  int ntri_file = 0;

  while (1) {

    // read a line and bcast length

    if (filereader) {
      if (fgets(line,MAXLINE,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
    }
    MPI_Bcast(&n,1,MPI_INT,0,filecomm);

    // if n = 0 then end-of-file so return with blank line

    if (n == 0) {
      line[0] = '\0';
      return;
    }

    // bcast line

    MPI_Bcast(line,n,MPI_CHAR,0,filecomm);

    // trim anything from '#' onward
    // if line is blank, continue

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    // search line for header keyword and set corresponding variable
    // else exit and line will store Section keyword

    if (strstr(line,"points")) {
      bigint bnpoint;
      sscanf(line,BIGINT_FORMAT,&bnpoint);
      if (bnpoint > MAXSMALLINT)
        error->one(FLERR,"Read surf npoint is too large");
      npoint_file = bnpoint;
    } else if (strstr(line,"lines")) {
      if (dim == 3)
        error->one(FLERR,"Surf file cannot contain lines for 3d simulation");
      bigint bnline;
      sscanf(line,BIGINT_FORMAT,&bnline);
      if (bnline > MAXSMALLINT) error->all(FLERR,"Read surf nline is too large");
      nsurf_file = nline_file = bnline;
    } else if (strstr(line,"triangles")) {
      if (dim == 2)
        error->one(FLERR,
                   "Surf file cannot contain triangles for 2d simulation");
      bigint bntri;
      sscanf(line,BIGINT_FORMAT,&bntri);
      if (bntri > MAXSMALLINT) error->one(FLERR,"Read surf ntri is too large");
      nsurf_file = ntri_file = bntri;

    } else break;
  }

  if (dim == 2 && nline_file == 0)
    error->one(FLERR,"Surf file does not contain lines");
  if (dim == 3 && ntri_file == 0)
    error->one(FLERR,"Surf file does not contain triangles");
}

/* ----------------------------------------------------------------------
   read all points from a single file
   store in local pts data struct until rest of file is read
------------------------------------------------------------------------- */

void ReadSurf::read_points()
{
  int i,m,n,nchunk;
  char *next,*buf;

  bigint nbytes = (bigint) npoint_file * sizeof(Point);
  pts = (Point *) memory->smalloc(nbytes,"readsurf:pts");

  // read and broadcast one CHUNK of points at a time

  int nread = 0;

  while (nread < npoint_file) {
    if (npoint_file-nread > CHUNK) nchunk = CHUNK;
    else nchunk = npoint_file - nread;
    if (filereader) {
      char *eof;
      m = 0;
      for (i = 0; i < nchunk; i++) {
        eof = fgets(&buffer[m],MAXLINE,fp);
        if (eof == NULL) error->one(FLERR,"Unexpected end of surf file");
        m += strlen(&buffer[m]);
      }
      m++;
    }
    MPI_Bcast(&m,1,MPI_INT,0,filecomm);
    MPI_Bcast(buffer,m,MPI_CHAR,0,filecomm);

    buf = buffer;
    next = strchr(buf,'\n');
    *next = '\0';
    int nwords = input->count_words(buf);
    *next = '\n';

    if (dim == 2 && nwords != 3)
      error->one(FLERR,"Incorrect point format in surf file");
    if (dim == 3 && nwords != 4)
      error->one(FLERR,"Incorrect point format in surf file");

    n = nread;

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
}

/* ----------------------------------------------------------------------
   read Lines section of file
------------------------------------------------------------------------- */

void ReadSurf::read_lines()
{
  int i,m,nchunk,type,p1,p2;
  surfint id;
  double x1[2],x2[2];
  char *next,*buf;

  // nwords_required = # of words per line

  int nwords_required;
  if (npoint_file) nwords_required = 3;
  else nwords_required = 5;
  if (typeflag) nwords_required++;
  nwords_required += nvalues_custom;
  double *custom = new double[nvalues_custom];

  // read and broadcast one CHUNK of lines at a time

  int nread = 0;

  while (nread < nsurf_file) {
    if (nsurf_file - nread > CHUNK) nchunk = CHUNK;
    else nchunk = nsurf_file - nread;
    if (filereader) {
      char *eof;
      m = 0;
      for (i = 0; i < nchunk; i++) {
        eof = fgets(&buffer[m],MAXLINE,fp);
        if (eof == NULL) error->one(FLERR,"Unexpected end of surf file");
        m += strlen(&buffer[m]);
      }
      m++;
    }
    MPI_Bcast(&m,1,MPI_INT,0,filecomm);
    MPI_Bcast(buffer,m,MPI_CHAR,0,filecomm);

    buf = buffer;
    next = strchr(buf,'\n');
    *next = '\0';
    int nwords = input->count_words(buf);
    *next = '\n';

    if (nwords != nwords_required)
      error->all(FLERR,"Incorrect line format in surf file");

    // if Points section in file, each read line has indices into it

    int i,ic,iv,icvalue;

    if (npoint_file) {
      for (i = 0; i < nchunk; i++) {
        next = strchr(buf,'\n');
        id = ATOSURFINT(strtok(buf," \t\n\r\f"));
        if (typeflag) type = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
        else type = 1;

        p1 = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
        p2 = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
        if (p1 < 1 || p1 > npoint_file || p2 < 1 || p2 > npoint_file || p1 == p2)
          error->all(FLERR,"Invalid point index in Lines section");

	if (ncustom) {
	  icvalue = 0;
	  for (ic = 0; ic < ncustom; ic++) {
	    if (type_custom[ic] == 0) {
	      if (size_custom[ic] == 0) {
                custom[icvalue++] =
                  input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
	      } else
		for (iv = 0; iv < size_custom[ic]; iv++)
		  custom[icvalue++] =
                    input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
	    } else {
              if (size_custom[ic] == 0)
		custom[icvalue++] =
                  input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
              else
		for (iv = 0; iv < size_custom[ic]; iv++)
		  custom[icvalue++] =
                    input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
            }
          }
	}

	if ((nread+i) % nprocs_file == me_file) {
	  add_line(id,type,pts[p1-1].x,pts[p2-1].x);
	  if (ncustom) add_custom(id,custom);
	  nsurf++;
	}

        buf = next + 1;
      }

    // if no Points section, each read line has point coords

    } else {
      for (int i = 0; i < nchunk; i++) {
        next = strchr(buf,'\n');
        id = ATOSURFINT(strtok(buf," \t\n\r\f"));
        if (typeflag) type = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
        else type = 1;

        x1[0] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
        x1[1] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
        x2[0] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
        x2[1] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
	
	if (ncustom) {
	  icvalue = 0;
	  for (ic = 0; ic < ncustom; ic++) {
	    if (type_custom[ic] == 0) {
	      if (size_custom[ic] == 0)
		custom[icvalue++] =
                  input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
	      else
		for (iv = 0; iv < size_custom[ic]; iv++)
		  custom[icvalue++] =
                    input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
	    } else {
	      if (size_custom[ic] == 0)
		custom[icvalue++] =
                  input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
	      else
		for (iv = 0; iv < size_custom[ic]; iv++)
		  custom[icvalue++] =
                    input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
            }
          }
	}

	if ((nread+i) % nprocs_file == me_file) {
	  add_line(id,type,x1,x2);
	  if (ncustom) add_custom(id,custom);
	  nsurf++;
	}

        buf = next + 1;
      }
    }

    nread += nchunk;
  }

  // clean up custom value storage

  delete [] custom;
}

/* ----------------------------------------------------------------------
   read Triangles section of file
------------------------------------------------------------------------- */

void ReadSurf::read_tris()
{
  int i,m,nchunk,type,p1,p2,p3;
  surfint id;
  double x1[3],x2[3],x3[3];
  char *next,*buf;

  // nwords_required = # of words per line

  int nwords_required;
  if (npoint_file) nwords_required = 4;
  else nwords_required = 10;
  if (typeflag) nwords_required++;
  nwords_required += nvalues_custom;
  double *custom = new double[nvalues_custom];

  // read and broadcast one CHUNK of triangles at a time

  int nread = 0;

  while (nread < nsurf_file) {
    if (nsurf_file - nread > CHUNK) nchunk = CHUNK;
    else nchunk = nsurf_file - nread;
    if (filereader) {
      char *eof;
      m = 0;
      for (i = 0; i < nchunk; i++) {
        eof = fgets(&buffer[m],MAXLINE,fp);
        if (eof == NULL) error->one(FLERR,"Unexpected end of surf file");
        m += strlen(&buffer[m]);
      }
      m++;
    }
    MPI_Bcast(&m,1,MPI_INT,0,filecomm);
    MPI_Bcast(buffer,m,MPI_CHAR,0,filecomm);

    buf = buffer;
    next = strchr(buf,'\n');
    *next = '\0';
    int nwords = input->count_words(buf);
    *next = '\n';

    if (nwords != nwords_required)
      error->all(FLERR,"Incorrect line format in surf file");

    // if Points section in file, each read line has indices into it

    int ic,iv,icvalue;

    if (npoint_file) {
      for (int i = 0; i < nchunk; i++) {
        next = strchr(buf,'\n');
        id = ATOSURFINT(strtok(buf," \t\n\r\f"));
        if (typeflag) type = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
        else type = 1;

        p1 = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
        p2 = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
        p3 = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
        if (p1 < 1 || p1 > npoint_file || p2 < 1 || p2 > npoint_file ||
            p3 < 1 || p3 > npoint_file || p1 == p2 || p2 == p3)
          error->all(FLERR,"Invalid point index in Triangles section");

	if (ncustom) {
	  icvalue = 0;
	  for (ic = 0; ic < ncustom; ic++) {
	    if (type_custom[ic] == 0) {
	      if (size_custom[ic] == 0)
		custom[icvalue++] =
                  input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
	      else
		for (iv = 0; iv < size_custom[ic]; iv++)
                  custom[icvalue++] =
                    input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
	    } else {
              if (size_custom[ic] == 0)
		custom[icvalue++] =
                  input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
	      else
		for (iv = 0; iv < size_custom[ic]; iv++)
		  custom[icvalue++] =
                    input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
            }
	  }
	}

	if ((nread+i) % nprocs_file == me_file) {
	  add_tri(id,type,pts[p1-1].x,pts[p2-1].x,pts[p3-1].x);
	  if (ncustom) add_custom(id,custom);
	  nsurf++;
	}
	
        buf = next + 1;
      }

    // if no Points section, each read line has point coords

    } else {
      for (int i = 0; i < nchunk; i++) {
        next = strchr(buf,'\n');
        id = ATOSURFINT(strtok(buf," \t\n\r\f"));
        if (typeflag) type = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
        else type = 1;

        x1[0] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
        x1[1] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
        x1[2] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
        x2[0] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
        x2[1] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
        x2[2] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
        x3[0] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
        x3[1] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
        x3[2] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));

	if (ncustom) {
	  icvalue = 0;
	  for (ic = 0; ic < ncustom; ic++) {
	    if (type_custom[ic] == 0) {
	      if (size_custom[ic] == 0)
		custom[icvalue++] = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
	      else
		for (iv = 0; iv < size_custom[ic]; iv++)
		  custom[icvalue++] = input->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
	    } else {
	      if (size_custom[ic] == 0)
		custom[icvalue++] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
	      else
		for (iv = 0; iv < size_custom[ic]; iv++)
		  custom[icvalue++] = input->numeric(FLERR,strtok(NULL," \t\n\r\f"));
            }
          }
	}

	if ((nread+i) % nprocs_file == me_file) {
	  add_tri(id,type,x1,x2,x3);
	  if (ncustom) add_custom(id,custom);
	  nsurf++;
	}
	
        buf = next + 1;
      }
    }

    nread += nchunk;
  }

  // clean up custom value storage

  delete [] custom;
}

/* ----------------------------------------------------------------------
   add a line to read-in lines stored by this proc
------------------------------------------------------------------------- */

void ReadSurf::add_line(surfint id, int itype, double *p1, double *p2)
{
  if (nsurf == maxsurf) {
    if ((bigint) maxsurf + DELTA > MAXSMALLINT)
      error->one(FLERR,"Read_surf add_line overflowed");
    maxsurf += DELTA;
    lines = (Surf::Line *)
      memory->srealloc(lines,maxsurf*sizeof(Surf::Line),"readsurf:lines");
    if (ncustom)
      memory->grow(cvalues,maxsurf,1+nvalues_custom,"readsurf:cvalues");
  }

  lines[nsurf].id = id;
  lines[nsurf].type = itype;
  lines[nsurf].mask = 1;
  lines[nsurf].isc = lines[nsurf].isr = -1;
  lines[nsurf].p1[0] = p1[0];
  lines[nsurf].p1[1] = p1[1];
  lines[nsurf].p1[2] = 0.0;
  lines[nsurf].p2[0] = p2[0];
  lines[nsurf].p2[1] = p2[1];
  lines[nsurf].p2[2] = 0.0;
  lines[nsurf].norm[0] = lines[nsurf].norm[1] = lines[nsurf].norm[2] = 0.0;
  lines[nsurf].transparent = 0;
}

/* ----------------------------------------------------------------------
   add a triangle to read-in triangles stored by this proc
------------------------------------------------------------------------- */

void ReadSurf::add_tri(surfint id, int itype, double *p1, double *p2, double *p3)
{
  if (nsurf == maxsurf) {
    if ((bigint) maxsurf + DELTA > MAXSMALLINT)
      error->one(FLERR,"Read_surf add_tri overflowed");
    maxsurf += DELTA;
    tris = (Surf::Tri *)
      memory->srealloc(tris,maxsurf*sizeof(Surf::Tri),"readsurf:tris");
    if (ncustom)
      memory->grow(cvalues,maxsurf,1+nvalues_custom,"readsurf:cvalues");
  }

  tris[nsurf].id = id;
  tris[nsurf].type = itype;
  tris[nsurf].mask = 1;
  tris[nsurf].isc = tris[nsurf].isr = -1;
  tris[nsurf].p1[0] = p1[0];
  tris[nsurf].p1[1] = p1[1];
  tris[nsurf].p1[2] = p1[2];
  tris[nsurf].p2[0] = p2[0];
  tris[nsurf].p2[1] = p2[1];
  tris[nsurf].p2[2] = p2[2];
  tris[nsurf].p3[0] = p3[0];
  tris[nsurf].p3[1] = p3[1];
  tris[nsurf].p3[2] = p3[2];
  tris[nsurf].norm[0] = tris[nsurf].norm[1] = tris[nsurf].norm[2] = 0.0;
  tris[nsurf].transparent = 0;
}

/* ----------------------------------------------------------------------
   add custom values for one line or one triangle to cvalues array
------------------------------------------------------------------------- */

void ReadSurf::add_custom(surfint id, double *custom)
{
  cvalues[nsurf][0] = ubuf(id).d;
  for (int ivalue = 0; ivalue < nvalues_custom; ivalue++)
    cvalues[nsurf][ivalue+1] = custom[ivalue];
}

// ---------------------------------------------------------------------
// methods to transform surface elements
// ---------------------------------------------------------------------

/* ----------------------------------------------------------------------
   parse additional keywords from start onward
   apply geometric transformation keywords immediately
   store keyword settngs for group/typeflag/transparent/partflag/file
     these are processed in constructor
------------------------------------------------------------------------- */

void ReadSurf::process_args(int start, int narg, char **arg)
{
  origin[0] = origin[1] = origin[2] = 0.0;
  grouparg = 0;
  typeadd = 0;
  partflag = NONE;
  transparent_flag = 0;
  filearg = 0;

  int iarg = start;
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

   // store these settings for later

    } else if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid read_surf command");
      grouparg = iarg+1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"typeadd") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid read_surf command");
      typeadd = input->inumeric(FLERR,arg[iarg+1]);
      if (typeadd < 0) error->all(FLERR,"Invalid read_surf command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"transparent") == 0) {
      transparent_flag = 1;
      iarg++;

    } else if (strcmp(arg[iarg],"particle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid read_surf command");
      if (strcmp(arg[iarg+1],"none") == 0) partflag = NONE;
      else if (strcmp(arg[iarg+1],"check") == 0) partflag = CHECK;
      else if (strcmp(arg[iarg+1],"keep") == 0) partflag = KEEP;
      else error->all(FLERR,"Invalid read_surf command");
      iarg += 2;

    // file must be last keyword, else WriteSurf will flag error

    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid read_surf command");
      filearg = iarg+1;
      iarg = narg;

    } else error->all(FLERR,"Invalid read_surf command");
  }

  // error test on particles

  if (particle->exist && partflag == NONE)
    error->all(FLERR,"Using read_surf particle none when particles exist");
}

/* ----------------------------------------------------------------------
   translate vertices by (dx,dy,dz)
   for 2d, dz will be 0.0
------------------------------------------------------------------------- */

void ReadSurf::translate(double dx, double dy, double dz)
{
  if (dim == 2) {
    for (int i = 0; i < nsurf; i++) {
      lines[i].p1[0] += dx;
      lines[i].p1[1] += dy;
      lines[i].p1[2] += dz;
      lines[i].p2[0] += dx;
      lines[i].p2[1] += dy;
      lines[i].p2[2] += dz;
    }

  } else if (dim == 3) {
    for (int i = 0; i < nsurf; i++) {
      tris[i].p1[0] += dx;
      tris[i].p1[1] += dy;
      tris[i].p1[2] += dz;
      tris[i].p2[0] += dx;
      tris[i].p2[1] += dy;
      tris[i].p2[2] += dz;
      tris[i].p3[0] += dx;
      tris[i].p3[1] += dy;
      tris[i].p3[2] += dz;
    }
  }
}

/* ----------------------------------------------------------------------
   scale vertices by (sx,sy,sz) around origin
   for 2d, do not reset x[2] to avoid epsilon change
------------------------------------------------------------------------- */

void ReadSurf::scale(double sx, double sy, double sz)
{
  if (dim == 2) {
    for (int i = 0; i < nsurf; i++) {
      lines[i].p1[0] = sx*(lines[i].p1[0]-origin[0]) + origin[0];
      lines[i].p1[1] = sy*(lines[i].p1[1]-origin[1]) + origin[1];
      lines[i].p2[0] = sx*(lines[i].p2[0]-origin[0]) + origin[0];
      lines[i].p2[1] = sy*(lines[i].p2[1]-origin[1]) + origin[1];
    }

  } else if (dim == 3) {
    for (int i = 0; i < nsurf; i++) {
      tris[i].p1[0] = sx*(tris[i].p1[0]-origin[0]) + origin[0];
      tris[i].p1[1] = sy*(tris[i].p1[1]-origin[1]) + origin[1];
      tris[i].p1[2] = sz*(tris[i].p1[2]-origin[2]) + origin[2];
      tris[i].p2[0] = sx*(tris[i].p2[0]-origin[0]) + origin[0];
      tris[i].p2[1] = sy*(tris[i].p2[1]-origin[1]) + origin[1];
      tris[i].p2[2] = sz*(tris[i].p2[2]-origin[2]) + origin[2];
      tris[i].p3[0] = sx*(tris[i].p3[0]-origin[0]) + origin[0];
      tris[i].p3[1] = sy*(tris[i].p3[1]-origin[1]) + origin[1];
      tris[i].p3[2] = sz*(tris[i].p3[2]-origin[2]) + origin[2];
    }
  }
}

/* ----------------------------------------------------------------------
   rotate vertices around origin
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

  if (dim == 2) {
    for (int i = 0; i < nsurf; i++) {
      d[0] = lines[i].p1[0] - origin[0];
      d[1] = lines[i].p1[1] - origin[1];
      d[2] = lines[i].p1[2] - origin[2];
      MathExtra::matvec(rotmat,d,dnew);
      lines[i].p1[0] = dnew[0] + origin[0];
      lines[i].p1[1] = dnew[1] + origin[1];

      d[0] = lines[i].p2[0] - origin[0];
      d[1] = lines[i].p2[1] - origin[1];
      d[2] = lines[i].p2[2] - origin[2];
      MathExtra::matvec(rotmat,d,dnew);
      lines[i].p2[0] = dnew[0] + origin[0];
      lines[i].p2[1] = dnew[1] + origin[1];
    }
  } else if (dim == 3) {
    for (int i = 0; i < nsurf; i++) {
      d[0] = tris[i].p1[0] - origin[0];
      d[1] = tris[i].p1[1] - origin[1];
      d[2] = tris[i].p1[2] - origin[2];
      MathExtra::matvec(rotmat,d,dnew);
      tris[i].p1[0] = dnew[0] + origin[0];
      tris[i].p1[1] = dnew[1] + origin[1];
      tris[i].p1[2] = dnew[2] + origin[2];

      d[0] = tris[i].p2[0] - origin[0];
      d[1] = tris[i].p2[1] - origin[1];
      d[2] = tris[i].p2[2] - origin[2];
      MathExtra::matvec(rotmat,d,dnew);
      tris[i].p2[0] = dnew[0] + origin[0];
      tris[i].p2[1] = dnew[1] + origin[1];
      tris[i].p2[2] = dnew[2] + origin[2];

      d[0] = tris[i].p3[0] - origin[0];
      d[1] = tris[i].p3[1] - origin[1];
      d[2] = tris[i].p3[2] - origin[2];
      MathExtra::matvec(rotmat,d,dnew);
      tris[i].p3[0] = dnew[0] + origin[0];
      tris[i].p3[1] = dnew[1] + origin[1];
      tris[i].p3[2] = dnew[2] + origin[2];
    }
  }
}

/* ----------------------------------------------------------------------
   invert vertex ordering within each line or tri
   this flips direction of surface normal
------------------------------------------------------------------------- */

void ReadSurf::invert()
{
  int tmp;
  double x[3];

  if (dim == 2) {
    for (int i = 0; i < nsurf; i++) {
      memcpy(x,lines[i].p1,3*sizeof(double));
      memcpy(lines[i].p1,lines[i].p2,3*sizeof(double));
      memcpy(lines[i].p2,x,3*sizeof(double));
    }

  } else if (dim == 3) {
    for (int i = 0; i < nsurf; i++) {
      memcpy(x,tris[i].p2,3*sizeof(double));
      memcpy(tris[i].p2,tris[i].p3,3*sizeof(double));
      memcpy(tris[i].p3,x,3*sizeof(double));
    }
  }
}

/* ----------------------------------------------------------------------
   clip all lines so fit inside simulation box
   may discard some lines completely
------------------------------------------------------------------------- */

void ReadSurf::clip2d()
{
  int i,dim,side,flag1,flag2;
  double value,param;
  double *p1,*p2,*inpt,*outpt;

  int *discard;
  memory->create(discard,nsurf,"readsurf:discard");
  for (i = 0; i < nsurf; i++) discard[i] = 0;
  int discardflag = 0;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  for (int iface = 0; iface < 4; iface++) {
    dim = iface / 2;
    side = iface % 2;
    if (side == 0) value = boxlo[dim];
    else value = boxhi[dim];

    // flag each point as (1,0,-1) = outside,on,inside clipping edge
    // keep line unchanged:
    //   if both pts are inside or on clipping edge
    //   at least one pt is inside clipping edge
    // discard line:
    //   if both pts are either outside or on clipping edge
    // straddle line:
    //   one pt is inside, one pt is outside
    //   replace outside pt with pt on the clipping edge

    for (i = 0; i < nsurf; i++) {
      if (discard[i]) continue;

      p1 = lines[i].p1;
      p2 = lines[i].p2;

      if (side == 0) {
        if (p1[dim] < value) flag1 = 1;
        else if (p1[dim] > value) flag1 = -1;
        else flag1 = 0;
        if (p2[dim] < value) flag2 = 1;
        else if (p2[dim] > value) flag2 = -1;
        else flag2 = 0;
      } else {
        if (p1[dim] > value) flag1 = 1;
        else if (p1[dim] < value) flag1 = -1;
        else flag1 = 0;
        if (p2[dim] > value) flag2 = 1;
        else if (p2[dim] < value) flag2 = -1;
        else flag2 = 0;
      }

      if (flag1 < 0 && flag2 <= 0) continue;
      if (flag1 <= 0 && flag2 < 0) continue;

      if (flag1 >= 0 && flag2 >= 0) {
        discard[i] = 1;
        discardflag = 1;
        continue;
      }

      if (flag1 < 0) {
        inpt = p1;
        outpt = p2;
      } else {
        inpt = p2;
        outpt = p1;
      }

      // reset outpt to intersection of line with clipping edge
      // last line insures outpt is exactly on clipping edge

      param = (value-inpt[dim]) / (outpt[dim]-inpt[dim]);
      outpt[0] = inpt[0] + param*(outpt[0]-inpt[0]);
      outpt[1] = inpt[1] + param*(outpt[1]-inpt[1]);
      outpt[dim] = value;
    }
  }

  // remove deleted lines

  int n = 0;
  for (i = 0; i < nsurf; i++) {
    if (!discard[i]) {
      if (n != i) {
        memcpy(&lines[n],&lines[i],sizeof(Surf::Line));
        if (ncustom)
          memcpy(cvalues[n],cvalues[i],(1+nvalues_custom)*sizeof(double));
      }
      n++;
    }
  }
  nsurf = n;

  // clean up

  memory->destroy(discard);

  // if discarded any surfs, renumber surf IDs across all procs

  int discardany;
  MPI_Allreduce(&discardflag,&discardany,1,MPI_INT,MPI_MAX,world);

  if (discardany) {
    bigint bnsurf = nsurf;
    bigint offset;
    MPI_Scan(&bnsurf,&offset,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    offset -= bnsurf;

    for (i = 0; i < nsurf; i++) {
      lines[i].id = static_cast<surfint> (offset + i + 1);
      if (ncustom) cvalues[i][0] = ubuf(lines[i].id).d;
    }

    MPI_Allreduce(&bnsurf,&nsurf_all,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  }

  // stats

  if (me == 0) {
    if (screen) fprintf(screen,"  clipped to " BIGINT_FORMAT " lines\n",nsurf_all);
    if (logfile) fprintf(logfile,"  clipped to " BIGINT_FORMAT " lines\n",nsurf_all);
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
  int i,dim,side,flag1,flag2,flag3,nin;
  double value,param;
  double x1[3],x2[3];
  double *p1,*p2,*p3,*in1,*in2,*out1,*out2;

  // discard flag for each surf
  // will be augmented in loop if tris are added

  int *discard;
  memory->create(discard,nsurf,"readsurf:discard");
  for (i = 0; i < nsurf; i++) discard[i] = 0;
  int maxdiscard = nsurf;
  int discardflag = 0;
  int addflag = 0;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  for (int iface = 0; iface < 6; iface++) {
    dim = iface / 2;
    side = iface % 2;
    if (side == 0) value = boxlo[dim];
    else value = boxhi[dim];

    // flag each point as (1,0,-1) = outside,on,inside clipping edge
    // keep tri unchanged:
    //   if all pts are inside or on clipping plane
    //   at least one pt is inside clipping plane
    // discard tri:
    //   if all pts are either outside or on clipping plane
    // straddle tri:
    //   at least one pt is inside, at least one pt is outside
    //   replace outside pts with pts on the clipping plane
    //   if 1 pt is inside, triangle remains
    //   if 2 pts are inside, trapezoid remains, convert to 2 tris
    //   latter requires adding a triangle to tris list
    // use nsurf_face for loop since nsurf may be incremented by add_tri()

    int nsurf_face = nsurf;

    for (i = 0; i < nsurf_face; i++) {
      if (discard[i]) continue;

      p1 = tris[i].p1;
      p2 = tris[i].p2;
      p3 = tris[i].p3;

      if (side == 0) {
        if (p1[dim] < value) flag1 = 1;
        else if (p1[dim] > value) flag1 = -1;
        else flag1 = 0;
        if (p2[dim] < value) flag2 = 1;
        else if (p2[dim] > value) flag2 = -1;
        else flag2 = 0;
        if (p3[dim] < value) flag3 = 1;
        else if (p3[dim] > value) flag3 = -1;
        else flag3 = 0;
      } else {
        if (p1[dim] > value) flag1 = 1;
        else if (p1[dim] < value) flag1 = -1;
        else flag1 = 0;
        if (p2[dim] > value) flag2 = 1;
        else if (p2[dim] < value) flag2 = -1;
        else flag2 = 0;
        if (p3[dim] > value) flag3 = 1;
        else if (p3[dim] < value) flag3 = -1;
        else flag3 = 0;
      }

      if (flag1 < 0 && flag2 <= 0 && flag3 <= 0) continue;
      if (flag1 <= 0 && flag2 < 0 && flag3 <= 0) continue;
      if (flag1 <= 0 && flag2 <= 0 && flag3 < 0) continue;

      if (flag1 >= 0 && flag2 >= 0 && flag3 >= 0) {
        discard[i] = 1;
        discardflag = 1;
        continue;
      }

      // nin = # of inside pts

      nin = 0;
      if (flag1 < 0) nin++;
      if (flag2 < 0) nin++;
      if (flag3 < 0) nin++;

      // one pt is inside
      // one or both of other 2 pts are outside
      // make a clip with each other pt that is outside
      // no new triangle is formed, existing tri is just changed

      if (nin == 1) {
        if (flag1 < 0) {
          in1 = p1;
          out1 = p2;
          out2 = p3;
        } else if (flag2 < 0) {
          in1 = p2;
          out1 = p3;
          out2 = p1;
        } else {
          in1 = p3;
          out1 = p1;
          out2 = p2;
        }
        if (out1[dim] != value) {
          param = (value-in1[dim]) / (out1[dim]-in1[dim]);
          out1[0] = in1[0] + param*(out1[0]-in1[0]);
          out1[1] = in1[1] + param*(out1[1]-in1[1]);
          out1[2] = in1[2] + param*(out1[2]-in1[2]);
          out1[dim] = value;
        }
        if (out2[dim] != value) {
          param = (value-in1[dim]) / (out2[dim]-in1[dim]);
          out2[0] = in1[0] + param*(out2[0]-in1[0]);
          out2[1] = in1[1] + param*(out2[1]-in1[1]);
          out2[2] = in1[2] + param*(out2[2]-in1[2]);
          out2[dim] = value;
        }
      }

      // two pts are inside, one pt is outside
      // set in1,in2,out1 so all 3 of these tris have same consistent normal
      // straddle tri = (in1,in2,out1)
      // modified tri = (in1,in2,x2)
      // added tri = (in1,x2,x1)
      // x1 = clip pt between in1 and out1
      // x2 = clip pt between in2 and out1
      // x1 and x2 will be computed exactly the same by a tri sharing the edge

      if (nin == 2) {
        if (flag1 > 0) {
          out1 = p1;
          in1 = p2;
          in2 = p3;
        } else if (flag2 > 0) {
          out1 = p2;
          in1 = p3;
          in2 = p1;
        } else {
          out1 = p3;
          in1 = p1;
          in2 = p2;
        }

        param = (value-in1[dim]) / (out1[dim]-in1[dim]);
        x1[0] = in1[0] + param*(out1[0]-in1[0]);
        x1[1] = in1[1] + param*(out1[1]-in1[1]);
        x1[2] = in1[2] + param*(out1[2]-in1[2]);
        x1[dim] = value;

        param = (value-in2[dim]) / (out1[dim]-in2[dim]);
        x2[0] = in2[0] + param*(out1[0]-in2[0]);
        x2[1] = in2[1] + param*(out1[1]-in2[1]);
        x2[2] = in2[2] + param*(out1[2]-in2[2]);
        x2[dim] = value;

        // reset one point in modified tri

        memcpy(out1,x2,3*sizeof(double));

        // add a new tri
        // use same ID as modified tri for now, will renumber below

        if (nsurf+1 >= maxdiscard) {
          maxdiscard += DELTA_DISCARD;
          memory->grow(discard,maxdiscard,"readsurf:discard");
        }

        // cannot use in1 in add_tri() because it points into tris,
        //   which may be realloced at start of add_tri()

        double in1_copy[3];
        in1_copy[0] = in1[0];
        in1_copy[1] = in1[1];
        in1_copy[2] = in1[2];

	add_tri(tris[i].id,tris[i].type,in1_copy,x2,x1);
	if (ncustom)
	  memcpy(cvalues[nsurf],cvalues[i],(1+nvalues_custom)*sizeof(double));
        discard[nsurf] = 0;
	nsurf++;
        addflag = 1;
      }
    }
  }

  // remove deleted tris

  int n = 0;
  for (i = 0; i < nsurf; i++) {
    if (!discard[i]) {
      if (n != i) {
        memcpy(&tris[n],&tris[i],sizeof(Surf::Tri));
        if (ncustom)
          memcpy(cvalues[n],cvalues[i],(1+nvalues_custom)*sizeof(double));
      }
      n++;
    }
  }
  nsurf = n;

  // clean up

  memory->destroy(discard);

  // if discarded or added any surfs, renumber surf IDs across all procs

  int changeflag = discardflag + addflag;
  int changeany;
  MPI_Allreduce(&changeflag,&changeany,1,MPI_INT,MPI_MAX,world);

  if (changeany) {
    bigint bnsurf = nsurf;
    bigint offset;
    MPI_Scan(&bnsurf,&offset,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    offset -= bnsurf;

    for (i = 0; i < nsurf; i++) {
      tris[i].id = static_cast<surfint> (offset + i + 1);
      if (ncustom) cvalues[i][0] = ubuf(tris[i].id).d;
    }

    MPI_Allreduce(&bnsurf,&nsurf_all,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  }

  // stats

  if (me == 0) {
    if (screen) fprintf(screen,"  clipped to " BIGINT_FORMAT " tris\n",nsurf_all);
    if (logfile) fprintf(logfile,"  clipped to " BIGINT_FORMAT " tris\n",nsurf_all);
  }
}

// ----------------------------------------------------------------------
// auxilliary methods
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   check that read-in surf IDs span bounds from 1 to N
   only checks if min/max surf IDs match bounds, not if truly consecutive
------------------------------------------------------------------------- */

void ReadSurf::check_bounds()
{
  surfint id;

  bigint smin = nsurf_all;
  bigint smax = 0;

  if (dim == 2) {
    for (int i = 0; i < nsurf; i++) {
      id = lines[i].id;
      smin = MIN(smin,id);
      smax = MAX(smax,id);
    }
  } else {
    for (int i = 0; i < nsurf; i++) {
      id = tris[i].id;
      smin = MIN(smin,id);
      smax = MAX(smax,id);
    }
  }

  bigint sminall,smaxall;
  MPI_Allreduce(&smin,&sminall,1,MPI_SPARTA_BIGINT,MPI_MIN,world);
  MPI_Allreduce(&smax,&smaxall,1,MPI_SPARTA_BIGINT,MPI_MAX,world);

  if (sminall != 1) {
    char str[128];
    sprintf(str,"Read_surf minimum surface ID is " BIGINT_FORMAT,sminall);
    error->all(FLERR,str);
  }

  if (smaxall != nsurf_all) {
    char str[128];
    sprintf(str,"Read_surf maximum surface ID is " BIGINT_FORMAT,smaxall);
    error->all(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   push all points to box boundary that are just inside
   1Dec20 added just outside logic
   delta = user-specified frac * (hi-lo)
   this avoids tiny clipped surf elements
------------------------------------------------------------------------- */

void ReadSurf::push_points_to_boundary(double frac)
{
  int i,j;
  double *x;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  double xdelta = frac * (boxhi[0]-boxlo[0]);
  double ydelta = frac * (boxhi[1]-boxlo[1]);
  double zdelta = frac * (boxhi[2]-boxlo[2]);

  if (dim == 2) {
    for (i = 0; i < nsurf; i++) {
      for (j = 0; j < 2; j++) {
        if (j == 0) x = lines[i].p1;
        else x = lines[i].p2;

        if (fabs(x[0]-boxlo[0]) < xdelta) x[0] = boxlo[0];
        if (fabs(x[0]-boxhi[0]) < xdelta) x[0] = boxhi[0];

        if (fabs(x[1]-boxlo[1]) < ydelta) x[1] = boxlo[1];
        if (fabs(x[1]-boxhi[1]) < ydelta) x[1] = boxhi[1];

        /*
        if (x[0] >= boxlo[0] && x[0] <= boxhi[0]) {
          if (x[0]-boxlo[0] < xdelta) x[0] = boxlo[0];
          else if (boxhi[0]-x[0] < xdelta) x[0] = boxhi[0];
        }
        if (x[1] >= boxlo[1] && x[1] <= boxhi[1]) {
          if (x[1]-boxlo[1] < ydelta) x[1] = boxlo[1];
          else if (boxhi[1]-x[1] < ydelta) x[1] = boxhi[1];
        }
        */
      }
    }

  } else if (dim == 3) {
    for (int i = 0; i < nsurf; i++) {
      for (j = 0; j < 3; j++) {
        if (j == 0) x = tris[i].p1;
        else if (j == 1) x = tris[i].p2;
        else x = tris[i].p3;

        if (fabs(x[0]-boxlo[0]) < xdelta) x[0] = boxlo[0];
        if (fabs(x[0]-boxhi[0]) < xdelta) x[0] = boxhi[0];

        if (fabs(x[1]-boxlo[1]) < ydelta) x[1] = boxlo[1];
        if (fabs(x[1]-boxhi[1]) < ydelta) x[1] = boxhi[1];

        if (fabs(x[2]-boxlo[2]) < zdelta) x[2] = boxlo[2];
        if (fabs(x[2]-boxhi[2]) < zdelta) x[2] = boxhi[2];

        /*
        if (x[0] >= boxlo[0] && x[0] <= boxhi[0]) {
          if (x[0]-boxlo[0] < xdelta) x[0] = boxlo[0];
          else if (boxhi[0]-x[0] < xdelta) x[0] = boxhi[0];
        }
        if (x[1] >= boxlo[1] && x[1] <= boxhi[1]) {
          if (x[1]-boxlo[1] < ydelta) x[1] = boxlo[1];
          else if (boxhi[1]-x[1] < ydelta) x[1] = boxhi[1];
        }
        if (x[2] >= boxlo[2] && x[2] <= boxhi[2]) {
          if (x[2]-boxlo[2] < zdelta) x[2] = boxlo[2];
          else if (boxhi[2]-x[2] < zdelta) x[2] = boxhi[2];
        }
        */
      }
    }
  }
}

/* ----------------------------------------------------------------------
   check norms of adjacent lines
   error if dot product of 2 norms is -1 -> infinitely thin surf
   warn if closer than EPSILON_NORM to -1
------------------------------------------------------------------------- */

void ReadSurf::check_neighbor_norm_2d()
{
  // NOTE: need to enable this for distributed
  // NOTE: and for all, now that surfs are stored with point coords
  // NOTE: may need rendezvous comm to do this
  // NOTE: could do this before transferring local lines to Surf::lines

  if (distributed || !distributed) return;

  /*
  int p1,p2;

  // count[I] = # of lines that vertex I is part of

  int *count;
  int **p2e;
  memory->create(count,npoint_file,"readsurf:count");
  memory->create(p2e,npoint_file,2,"readsurf:count");
  for (int i = 0; i < npoint_file; i++) count[i] = 0;

  for (int i = 0; i < nsurf_file; i++) {
    //p1 = lines[i].p1;
    //p2 = lines[i].p2;
    p2e[p1][count[p1]++] = i;
    p2e[p2][count[p2]++] = i;
  }

  // check that norms of adjacent lines are not in opposite directions
  // norms are stored in Surf::lines, at end of orignal nsurf_old surfs

  Surf::Line *surflines = surf->lines;

  double dot;
  double *norm1,*norm2;

  int nerror = 0;
  int nwarn = 0;
  for (int i = 0; i < npoint_file; i++) {
    if (count[i] == 1) continue;
    norm1 = surflines[p2e[i][0] + nsurf_old].norm;
    norm2 = surflines[p2e[i][1] + nsurf_old].norm;
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
  */
}

/* ----------------------------------------------------------------------
   check norms of new adjacent triangles
   error if dot product of 2 norms is -1 -> infinitely thin surf
   warn if closer than EPSILON_NORM to -1
------------------------------------------------------------------------- */

void ReadSurf::check_neighbor_norm_3d()
{
  // NOTE: need to enable this for distributed
  // NOTE: and for all, now that surfs are stored with point coords
  // NOTE: may need rendezvous comm to do this
  // NOTE: could do this before transferring local tris to Surf::tris

  if (distributed || !distributed) return;

  /*
  int ntri_file;

  // hash directed edges of all triangles
  // key = directed edge, value = triangle it is part of
  // NOTE: could prealloc hash to correct size here

  MyHash hash;
  MyIterator it;

  // insert each edge into hash with triangle index as value

  bigint p1,p2,p3,key;

  for (int i = 0; i < ntri_file; i++) {
    //p1 = tris[i].p1;
    //p2 = tris[i].p2;
    //p3 = tris[i].p3;
    key = (p1 << 32) | p2;
    hash[key] = i;
    key = (p2 << 32) | p3;
    hash[key] = i;
    key = (p3 << 32) | p1;
    hash[key] = i;
  }

  // check that norms of adjacent triangles are not in opposite directions
  // norms are stored in Surf::tris, at end of orignal nsurf_old surfs

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
    norm1 = surftris[it->second + nsurf_old].norm;
    norm2 = surftris[hash[key] + nsurf_old].norm;
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
  */
}

/* ----------------------------------------------------------------------
   reading proc opens data file
   test if gzipped
   sets fp = file pointer
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
   infile contains a "*"
   search for all files which match the infile pattern
   replace "*" with latest timestep value to create outfile name
   search dir referenced by initial pathname of file
   if infile also contains "%", use "base" when searching directory
   only called by proc 0
------------------------------------------------------------------------- */

void ReadSurf::file_search(char *infile, char *outfile)
{
  char *ptr;

  // separate infile into dir + filename

  char *dirname = new char[strlen(infile) + 1];
  char *filename = new char[strlen(infile) + 1];

  if (strchr(infile,'/')) {
    ptr = strrchr(infile,'/');
    *ptr = '\0';
    strcpy(dirname,infile);
    strcpy(filename,ptr+1);
    *ptr = '/';
  } else {
    strcpy(dirname,"./");
    strcpy(filename,infile);
  }

  // if filename contains "%" replace "%" with "base"

  char *pattern = new char[strlen(filename) + 16];

  if ((ptr = strchr(filename,'%'))) {
    *ptr = '\0';
    sprintf(pattern,"%s%s%s",filename,"base",ptr+1);
    *ptr = '%';
  } else strcpy(pattern,filename);

  // scan all files in directory, searching for files that match pattern
  // maxnum = largest int that matches "*"

  size_t n = strlen(pattern) + 16;
  char *begin = new char[n];
  char *middle = new char[n];
  char *end = new char[n];

  ptr = strchr(pattern,'*');
  *ptr = '\0';
  strcpy(begin,pattern);
  strcpy(end,ptr+1);
  int nbegin = strlen(begin);
  bigint maxnum = -1;

  struct dirent *ep;
  DIR *dp = opendir(dirname);
  if (dp == NULL)
    error->one(FLERR,"Cannot open dir to search for restart file");
  while ((ep = readdir(dp))) {
    if (strstr(ep->d_name,begin) != ep->d_name) continue;
    if ((ptr = strstr(&ep->d_name[nbegin],end)) == NULL) continue;
    if (strlen(end) == 0) ptr = ep->d_name + strlen(ep->d_name);
    *ptr = '\0';
    if (strlen(&ep->d_name[nbegin]) < n) {
      strcpy(middle,&ep->d_name[nbegin]);
      if (ATOBIGINT(middle) > maxnum) maxnum = ATOBIGINT(middle);
    }
  }
  closedir(dp);
  if (maxnum < 0) error->one(FLERR,"Found no restart file matching pattern");

  // create outfile with maxint substituted for "*"
  // use original infile, not pattern, since need to retain "%" in filename

  ptr = strchr(infile,'*');
  *ptr = '\0';
  sprintf(outfile,"%s" BIGINT_FORMAT "%s",infile,maxnum,ptr+1);
  *ptr = '*';

  // clean up

  delete [] dirname;
  delete [] filename;
  delete [] pattern;
  delete [] begin;
  delete [] middle;
  delete [] end;
}

/* ----------------------------------------------------------------------
   grab next keyword
   if first = 1, line holds non-blank line that ended header
   else read lines until one is non-blank
   keyword is all text on line w/out leading & trailing white space
   read one additional line (assumed blank)
   if any read hits EOF, set keyword to empty string
------------------------------------------------------------------------- */

void ReadSurf::parse_keyword(int first)
{
  int eof = 0;

  // filereader reads upto non-blank line plus 1 following line
  // eof is set to 1 if any read hits end-of-file

  if (filereader) {
    if (!first) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
    }
    while (eof == 0 && strspn(line," \t\n\r") == strlen(line)) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
    }
    if (fgets(buffer,MAXLINE,fp) == NULL) eof = 1;
  }

  // if eof, set keyword empty and return

  MPI_Bcast(&eof,1,MPI_INT,0,filecomm);
  if (eof) {
    keyword[0] = '\0';
    return;
  }

  // bcast keyword line to all procs

  int n;
  if (filereader) n = strlen(line) + 1;
  MPI_Bcast(&n,1,MPI_INT,0,filecomm);
  MPI_Bcast(line,n,MPI_CHAR,0,filecomm);

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
