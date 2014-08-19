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

#include "stdlib.h"
#include "string.h"
#include "read_grid.h"
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

#define MAXLINE 256
#define CHUNK 1024

/* ---------------------------------------------------------------------- */

ReadGrid::ReadGrid(SPARTA *sparta) : Pointers(sparta)
{
  line = new char[MAXLINE];
  keyword = new char[MAXLINE];
  buffer = new char[CHUNK*MAXLINE];
}

/* ---------------------------------------------------------------------- */

ReadGrid::~ReadGrid()
{
  delete [] line;
  delete [] keyword;
  delete [] buffer;
}

/* ----------------------------------------------------------------------
   called as read_grid command in input script
------------------------------------------------------------------------- */

void ReadGrid::command(int narg, char **arg)
{
  if (!domain->box_exist) 
    error->all(FLERR,"Cannot read grid before simulation box is defined");
  if (grid->exist)
    error->all(FLERR,"Cannot read grid when grid is already defined");

  grid->exist = 1;

  if (narg != 1) error->all(FLERR,"Illegal read_grid command");

  read(arg[0],0);
}

/* ----------------------------------------------------------------------
   called from command() and directly from AdaptGrid
------------------------------------------------------------------------- */

void ReadGrid::read(char *filename, int external)
{
  // read file, create parent cells and then child cells

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  me = comm->me;
  if (me == 0) {
    if (!external && screen) 
      fprintf(screen,"Reading grid file ... %s\n",filename);
    open(filename);
  }

  header();

  // clear Grid::hash before re-populating it

#ifdef SPARTA_MAP
  std::map<cellint,int> *hash = grid->hash;
#else
  std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

  hash->clear();

  // read Parents section
  // create one parent cell per line

  parse_keyword(1);
  if (strcmp(keyword,"Parents") != 0)
    error->all(FLERR,
	       "Read_grid did not find parents section of grid file");
  read_parents();

  // close file

  if (me == 0) {
    if (compressed) pclose(fp);
    else fclose(fp);
  }

  // induce child cells from parent cells

  create_children();

  // clear Grid::hash since done using it

  hash->clear();
  grid->hashfilled = 0;

  // invoke grid methods to complete grid setup

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->find_neighbors();
  grid->check_uniform();
  comm->reset_neighbors();

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // stats

  if (external) return;

  double time_total = time3-time1;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  child cells = " BIGINT_FORMAT "\n",grid->ncell);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  read/ghost percent = %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total);
    }

    if (logfile) {
      fprintf(logfile,"  child cells = " BIGINT_FORMAT "\n",grid->ncell);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  read/ghost percent = %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   create one parent cell per line of Parents section of grid file
------------------------------------------------------------------------- */

void ReadGrid::create_parents(int n, char *buf)
{
  int j,m;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = count_words(buf);
  *next = '\n';

  if (nwords != 5)
    error->all(FLERR,"Incorrect format of parent cell in grid file");

  char **values = new char*[nwords];

  // loop over lines of parents
  // tokenize the line into values
  // create one parent cell for each line

  cellint id,idparent,ichild;
  int iparent,nx,ny,nz;
  char pstr[32],cstr[32];
  double lo[3],hi[3];
  Grid::ParentCell *p;

  // use Grid hash to store key = parent ID, value = index in pcells
  // NOTE: could prealloc hash to size nparents here

#ifdef SPARTA_MAP
  std::map<cellint,int> *hash = grid->hash;
#else
  std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

  Grid::ParentCell *pcells = grid->pcells;
  int dimension = domain->dimension;

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    for (j = 1; j < nwords; j++)
      values[j] = strtok(NULL," \t\n\r\f");

    // convert values[1] into id
    // create parent cell from id
    // error if parent cell already exists or its parent does not
    // NOTE: add more error checks on values[1]

    id = grid->id_str2num(values[1]);
    if (id < 0) error->all(FLERR,"Invalid cell ID in grid file");

    if (id) {
      if (hash->find(id) != hash->end()) 
        error->all(FLERR,"Duplicate cell ID in grid file");
      grid->id_pc_split(values[1],pstr,cstr);
      idparent = grid->id_str2num(pstr);
      ichild = ATOCELLINT(cstr);
      iparent = (*hash)[idparent];
      if (iparent < 0) 
        error->all(FLERR,"Parent cell's parent does not exist in grid file");
    } else iparent = -1;

    nx = atoi(values[2]);
    ny = atoi(values[3]);
    nz = atoi(values[4]);
    if (nx <= 0 || ny <= 0 || nz <= 0)
      error->all(FLERR,"Invalid Nx,Ny,Nz values in grid file");
    if (dimension == 2 && nz != 1) 
      error->all(FLERR,"Nz value in read_grid file must be 1 "
                 "for a 2d simulation");

    if (id) grid->id_child_lohi(iparent,ichild,lo,hi);
    else {
      lo[0] = domain->boxlo[0]; 
      lo[1] = domain->boxlo[1];
      lo[2] = domain->boxlo[2];
      hi[0] = domain->boxhi[0]; 
      hi[1] = domain->boxhi[1];
      hi[2] = domain->boxhi[2];
    } 

    grid->add_parent_cell(id,iparent,nx,ny,nz,lo,hi);
    (*hash)[id] = grid->nparent - 1;

    buf = next + 1;
  }

  delete [] values;
}

/* ----------------------------------------------------------------------
   create child cells from parent cells whose children are not parents
------------------------------------------------------------------------- */

void ReadGrid::create_children()
{
  Grid::ParentCell *pcells = grid->pcells;
  int nparent = grid->nparent;

#ifdef SPARTA_MAP
  std::map<cellint,int> *hash = grid->hash;
#else
  std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

  int nprocs = comm->nprocs;
  bigint count = 0;

  // loop over parent cells and its 3d array of child cells
  // if a child does not exist as a parent, create it as a child cell
  // assign child cells to procs in round-robin fashion via count

  // NOTE: change this to assign child to proc based on mod of cell ID
  //       could do this only when adapt grid invokes it

  int ix,iy,iz,nx,ny,nz,nbits;
  cellint m,idparent,idchild;
  double lo[3],hi[3];
  Grid::ParentCell *p;

  for (int iparent = 0; iparent < nparent; iparent++) {
    p = &pcells[iparent];
    idparent = p->id;
    nbits = p->nbits;
    nx = p->nx;
    ny = p->ny;
    nz = p->nz;

    m = 0;
    for (iz = 0; iz < nz; iz++)
      for (iy = 0; iy < ny; iy++)
        for (ix = 0; ix < nx; ix++) {
          m++;
          idchild = idparent | (m << nbits);
          if (hash->find(idchild) == hash->end()) {
            if (count % nprocs == me) {
              grid->id_child_lohi(iparent,m,lo,hi);
              grid->add_child_cell(idchild,iparent,lo,hi);
            }
            count++;
          }
        }
  }
}

/* ----------------------------------------------------------------------
   proc 0 opens grid file
   test if gzipped
------------------------------------------------------------------------- */

void ReadGrid::open(char *file)
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
   read free-format header of grid file
   1st line and blank lines are skipped
   non-blank lines are checked for header keywords and leading value is read
   header ends with non-blank line containing no header keyword (or EOF)
   return line with non-blank line (or empty line if EOF)
------------------------------------------------------------------------- */

void ReadGrid::header()
{
  int n;
  char *ptr;

  // skip 1st line of file

  if (me == 0) {
    char *eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of grid file");
  }

  nparents = 0;

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

    if (strstr(line,"parents")) sscanf(line,"%d",&nparents);
    else break;
  }

  if (nparents == 0) error->all(FLERR,"Grid file does not contain parents");
}

/* ----------------------------------------------------------------------
   grab next keyword
   read lines until one is non-blank
   keyword is all text on line w/out leading & trailing white space
   read one additional line (assumed blank)
   if any read hits EOF, set keyword to empty
   if first = 1, line variable holds non-blank line that ended header
------------------------------------------------------------------------- */

void ReadGrid::parse_keyword(int first)
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

int ReadGrid::count_words(char *line)
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
   read/store all parents
------------------------------------------------------------------------- */

void ReadGrid::read_parents()
{
  int i,m,nchunk;
  char *next,*buf;

  // read and broadcast one CHUNK of lines at a time

  int n = 0;
  int nread = 0;
  
  while (nread < nparents) {
    if (nparents-nread > CHUNK) nchunk = CHUNK;
    else nchunk = nparents-nread;
    if (me == 0) {
      char *eof;
      m = 0;
      for (i = 0; i < nchunk; i++) {
	eof = fgets(&buffer[m],MAXLINE,fp);
	if (eof == NULL) error->one(FLERR,"Unexpected end of grid file");
	m += strlen(&buffer[m]);
      }
      if (buffer[m-1] != '\n') strcpy(&buffer[m++],"\n");
      m++;
    }
    MPI_Bcast(&m,1,MPI_INT,0,world);
    MPI_Bcast(buffer,m,MPI_CHAR,0,world);

    create_parents(nchunk,buffer);
    nread += nchunk;
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  %d parent cells\n",nparents);
    if (logfile) fprintf(logfile,"  %d parent cells\n",nparents);
  }
}
