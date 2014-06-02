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

#include "spatype.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "dirent.h"
#include "read_restart.h"
#include "universe.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "grid.h"
#include "surf.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

// same as write_restart.cpp

#define MAGIC_STRING "SpartA RestartT"
#define ENDIAN 0x0001
#define ENDIANSWAP 0x1000
#define VERSION_NUMERIC 0

enum{VERSION,SMALLINT,CELLINT,BIGINT,
     UNITS,NTIMESTEP,NPROCS,
     FNUM,NRHO,VSTREAM,TEMP_THERMAL,GRAVITY,SURFMAX,GRIDCUT,
     COMM_SORT,COMM_STYLE,
     DIMENSION,AXISYMMETRIC,BOXLO,BOXHI,BFLAG,
     NPARTICLE,NUNSPLIT,NSPLIT,NSUB,NPOINT,NSURF,
     PARTICLE,GRID,SURF,
     MULTIPROC,PROCSPERFILE,PERPROC};

/* ---------------------------------------------------------------------- */

ReadRestart::ReadRestart(SPARTA *spa) : Pointers(spa) {}

/* ---------------------------------------------------------------------- */

void ReadRestart::command(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal read_restart command");

  if (domain->box_exist)
    error->all(FLERR,"Cannot read_restart after simulation box is defined");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // if filename contains "*", search dir for latest restart file

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

  // open single restart file or base file for multiproc case

  if (me == 0) {
    if (screen) fprintf(screen,"Reading restart file ...\n");
    char *hfile;
    if (multiproc) {
      hfile = new char[strlen(file) + 16];
      char *ptr = strchr(file,'%');
      *ptr = '\0';
      sprintf(hfile,"%s%s%s",file,"base",ptr+1);
      *ptr = '%';
    } else hfile = file;
    fp = fopen(hfile,"rb");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open restart file %s",hfile);
      error->one(FLERR,str);
    }
    if (multiproc) delete [] hfile;
  }

  // read magic string, endian flag, numeric version

  magic_string();
  endian();
  int incompatible = version_numeric();

  // read header info which creates simulation box
  // also defines species, grid, surfs

  header(incompatible);

  box_params();
  domain->box_exist = 1;
  particle_params();
  grid_params();
  grid->exist = 1;
  surf->exist = surf_params();

  // read file layout info

  file_layout();

  // close header file if in multiproc mode

  if (multiproc && me == 0) fclose(fp);

  // add parent cells to Grid::hash

#ifdef SPARTA_MAP
  std::map<cellint,int> *hash = grid->hash;
#else
  std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

  Grid::ParentCell *pcells = grid->pcells;
  int nparent = grid->nparent;

  hash->clear();
  for (int icell = 0; icell < nparent; icell++)
    (*hash)[pcells[icell].id] = -(icell+1);

  // read per-proc info, grid cells and particles

  int m,n,flag,value,tmp;
  long filepos;
  MPI_Status status;
  MPI_Request request;

  int maxbuf = 0;
  char *buf = NULL;

  // input of single native file
  // same proc count in file and current simulation
  // each proc will own exactly what it owned in previous run
  // proc 0 reads a chunk and sends it to owning proc
  // except skips its own chunk, then reads it at end
  // each proc:
  //   creates its grid cells from cell IDs
  //   assigns all particles to its cells

  if (multiproc == 0 && nprocs_file == nprocs) {

    if (me == 0) {
      for (int iproc = 0; iproc < nprocs_file; iproc++) {
        fread(&value,sizeof(int),1,fp);
        if (value != PERPROC)
          error->one(FLERR,"Invalid flag in peratom section of restart file");

        if (iproc == 0) filepos = ftell(fp);

        fread(&n,sizeof(int),1,fp);
        if (n > maxbuf) {
          maxbuf = n;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }

        if (iproc > 0) {
          fread(buf,sizeof(char),n,fp);
          MPI_Send(&n,1,MPI_INT,iproc,0,world);
          MPI_Recv(&tmp,0,MPI_INT,iproc,0,world,&status);
          MPI_Send(buf,n,MPI_CHAR,iproc,0,world);
        } else fseek(fp,filepos+sizeof(int)+n,SEEK_SET);
      }

      // rewind and read my chunk

      fseek(fp,filepos,SEEK_SET);
      fread(&n,sizeof(int),1,fp);
      fread(buf,sizeof(char),n,fp);

      fclose(fp);

    } else {
      MPI_Recv(&n,1,MPI_INT,0,0,world,&status);
      if (n > maxbuf) {
        maxbuf = n;
        memory->destroy(buf);
        memory->create(buf,maxbuf,"read_restart:buf");
      }
      MPI_Irecv(buf,n,MPI_CHAR,0,0,world,&request);
      MPI_Send(&tmp,0,MPI_INT,0,0,world);
      MPI_Wait(&request,&status);
    }

    n = grid->unpack_restart(buf);
    create_child_cells(0);
    n += particle->unpack_restart(&buf[n]);
    assign_particles(0);
  }

  // input of single native file
  // different proc count in file and current simulation
  // proc 0 reads a chunk and bcasts it to all other procs
  // each proc:
  //   creates 1/P fraction of grid cells from cell IDs
  //   assigns all particles to cells it owns

  else if (multiproc == 0) {

    for (int iproc = 0; iproc < nprocs_file; iproc++) {
      value = read_int();
      if (value != PERPROC)
        error->one(FLERR,"Invalid flag in peratom section of restart file");

      n = read_int();
      if (n > maxbuf) {
        maxbuf = n;
        memory->destroy(buf);
        memory->create(buf,maxbuf,"read_restart:buf");
      }

      read_char_vec(n,buf);

      n = grid->unpack_restart(buf);
      create_child_cells(1);
      n += particle->unpack_restart(&buf[n]);
      assign_particles(1);
    }
  }

  // input of multiple native files with procs <= files
  // # of files = multiproc_file
  // each proc reads a subset of files, striding by nprocs
  // each proc keeps all cells/particles in all perproc chunks in its files

  else if (nprocs <= multiproc_file) {

    char *procfile = new char[strlen(file) + 16];
    char *ptr = strchr(file,'%');

    for (int iproc = me; iproc < multiproc_file; iproc += nprocs) {
      *ptr = '\0';
      sprintf(procfile,"%s%d%s",file,iproc,ptr+1);
      *ptr = '%';
      fp = fopen(procfile,"rb");
      if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open restart file %s",procfile);
        error->one(FLERR,str);
      }

      fread(&flag,sizeof(int),1,fp);
      if (flag != PROCSPERFILE) 
        error->one(FLERR,"Invalid flag in peratom section of restart file");
      int procsperfile;
      fread(&procsperfile,sizeof(int),1,fp);

      for (int i = 0; i < procsperfile; i++) {
        fread(&flag,sizeof(int),1,fp);
        if (flag != PERPROC) 
          error->one(FLERR,"Invalid flag in peratom section of restart file");
        
        fread(&n,sizeof(int),1,fp);
        if (n > maxbuf) {
          maxbuf = n;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }
        fread(buf,sizeof(char),n,fp);

        n = grid->unpack_restart(buf);
        create_child_cells(0);
        n += particle->unpack_restart(&buf[n]);
        assign_particles(0);
      }

      fclose(fp);
    }

    delete [] procfile;
  }

  // input of multiple native files with procs > files
  // # of files = multiproc_file
  // cluster procs based on # of files
  // 1st proc in each cluster reads per-proc chunks from file
  // sends chunks round-robin to other procs in its cluster
  // each proc keeps all cells/particles in its perproc chunks in file

  else {

    // nclusterprocs = # of procs in my cluster that read from one file
    // filewriter = 1 if this proc reads file, else 0
    // fileproc = ID of proc in my cluster who reads from file
    // clustercomm = MPI communicator within my cluster of procs

    int nfile = multiproc_file;
    int icluster = static_cast<int> ((bigint) me * nfile/nprocs);
    int fileproc = static_cast<int> ((bigint) icluster * nprocs/nfile);
    int fcluster = static_cast<int> ((bigint) fileproc * nfile/nprocs);
    if (fcluster < icluster) fileproc++;
    int fileprocnext = 
      static_cast<int> ((bigint) (icluster+1) * nprocs/nfile);
    fcluster = static_cast<int> ((bigint) fileprocnext * nfile/nprocs);
    if (fcluster < icluster+1) fileprocnext++;
    int nclusterprocs = fileprocnext - fileproc;
    int filereader = 0;
    if (me == fileproc) filereader = 1;
    MPI_Comm clustercomm;
    MPI_Comm_split(world,icluster,0,&clustercomm);

    if (filereader) {
      char *procfile = new char[strlen(file) + 16];
      char *ptr = strchr(file,'%');
      *ptr = '\0';
      sprintf(procfile,"%s%d%s",file,icluster,ptr+1);
      *ptr = '%';
      fp = fopen(procfile,"rb");
      if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open restart file %s",procfile);
        error->one(FLERR,str);
      }
      delete [] procfile;
    }

    int flag,procsperfile;

    if (filereader) {
      fread(&flag,sizeof(int),1,fp);
      if (flag != PROCSPERFILE) 
        error->one(FLERR,"Invalid flag in peratom section of restart file");
      fread(&procsperfile,sizeof(int),1,fp);
    }
    MPI_Bcast(&procsperfile,1,MPI_INT,0,clustercomm);

    int tmp,iproc;
    MPI_Status status;
    MPI_Request request;

    for (int i = 0; i < procsperfile; i++) {
      if (filereader) {
        fread(&flag,sizeof(int),1,fp);
        if (flag != PERPROC) 
          error->one(FLERR,"Invalid flag in peratom section of restart file");

        fread(&n,sizeof(int),1,fp);
        if (n > maxbuf) {
          maxbuf = n;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }
        fread(buf,sizeof(char),n,fp);

        if (i % nclusterprocs) {
          iproc = me + (i % nclusterprocs);
          MPI_Send(&n,1,MPI_INT,iproc,0,world);
          MPI_Recv(&tmp,0,MPI_INT,iproc,0,world,&status);
          MPI_Rsend(buf,n,MPI_CHAR,iproc,0,world);
        }

      } else if (i % nclusterprocs == me - fileproc) {
        MPI_Recv(&n,1,MPI_INT,fileproc,0,world,&status);
        if (n > maxbuf) {
          maxbuf = n;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }
        MPI_Irecv(buf,n,MPI_CHAR,fileproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,fileproc,0,world);
        MPI_Wait(&request,&status);
      }

      if (i % nclusterprocs == me - fileproc) {
        n = grid->unpack_restart(buf);
        create_child_cells(0);
        n += particle->unpack_restart(&buf[n]);
        assign_particles(0);
      }
    }

    if (filereader) fclose(fp);
    MPI_Comm_free(&clustercomm);
  }

  // clean-up memory

  delete [] file;
  memory->destroy(buf);

  // clear Grid::hash since done using it

  hash->clear();
  grid->hashfilled = 0;

  // grid is no longer clumped unless reading on same # of procs
  // clumped decomposition is maintained (assuiming original file had it)
  //   for all reading methods above where nprocs_file = current nprocs

  if (nprocs_file != nprocs) grid->clumped = 0;

  // invoke surf and grid methods to complete grid setup

  if (surf->exist) {
    surf->setup_surf();
    grid->surf2grid(1);
  }

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->find_neighbors();
  grid->check_uniform();
  comm->reset_neighbors();

  if (surf->exist) {
    grid->set_inout();
    grid->type_check();
  }

  // check that all grid cells and particles were assigned to procs
  // print stats on grid cells, particles, surfs

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " grid cells\n",
                        grid->ncell);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " grid cells\n",
                        grid->ncell);
  }

  if (grid->nunsplit != nunsplit_file) 
    error->all(FLERR,"Did not assign all restart unsplit grid cells correctly");
  if (grid->nsplit != nsplit_file) 
    error->all(FLERR,"Did not assign all restart split grid cells correctly");
  if (grid->nsub != nsub_file) 
    error->all(FLERR,"Did not assign all restart sub grid cells correctly");

  bigint btmp = particle->nlocal;
  MPI_Allreduce(&btmp,&particle->nglobal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " particles\n",
                        particle->nglobal);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " particles\n",
                         particle->nglobal);
  }

  if (particle->nglobal != nparticle_file) 
    error->all(FLERR,"Did not assign all restart particles correctly");

  if (me == 0 && surf->exist) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " surf points\n",
                        surf->npoint);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " surf points\n",
                         surf->npoint);
    if (domain->dimension == 2) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " surf lines\n",
                          surf->nline);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " surf lines\n",
                           surf->nline);
    } else {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " surf triangles\n",
                          surf->ntri);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " surf triangles\n",
                           surf->ntri);
    }
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

void ReadRestart::file_search(char *infile, char *outfile)
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

  int n = strlen(pattern) + 16;
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
   read header of restart file
------------------------------------------------------------------------- */

void ReadRestart::header(int incompatible)
{
  // read flags and fields until flag = -1

  int flag = read_int();
  while (flag >= 0) {

    // check restart file version, warn if different

    if (flag == VERSION) {
      char *version = read_string();
      if (me == 0) {
        if (screen) fprintf(screen,"  restart file = %s, SPARTA = %s\n",
                            version,universe->version);
      }
      if (incompatible) 
        error->all(FLERR,"Restart file incompatible with current version");
      delete [] version;

    // check spatype.h sizes, error if different

    } else if (flag == SMALLINT) {
      int size = read_int();
      if (size != sizeof(smallint))
        error->all(FLERR,"Smallint setting in spatype.h is not compatible");
    } else if (flag == CELLINT) {
      int size = read_int();
      if (size != sizeof(cellint))
        error->all(FLERR,"Cellint setting in spatype.h is not compatible");
    } else if (flag == BIGINT) {
      int size = read_int();
      if (size != sizeof(bigint))
        error->all(FLERR,"Bigint setting in spatype.h is not compatible");

    // reset unit_style only if different
    // so that timestep is not changed

    } else if (flag == UNITS) {
      char *style = read_string();
      if (strcmp(style,update->unit_style) != 0) update->set_units(style);
      delete [] style;

    } else if (flag == NTIMESTEP) {
      update->ntimestep = read_bigint();

    // read nprocs from restart file, warn if different

    } else if (flag == NPROCS) {
      nprocs_file = read_int();
      if (nprocs_file != comm->nprocs && me == 0)
        error->warning(FLERR,"Restart file used different # of processors");

    // global settings

    } else if (flag == FNUM) {
      update->fnum = read_double();
    } else if (flag == NRHO) {
      update->nrho = read_double();
    } else if (flag == VSTREAM) {
      read_int();
      read_double_vec(3,update->vstream);
    } else if (flag == TEMP_THERMAL) {
      update->temp_thermal = read_double();
    } else if (flag == GRAVITY) {
      read_int();
      read_double_vec(3,update->gravity);
    } else if (flag == SURFMAX) {
      grid->maxsurfpercell = read_int();
    } else if (flag == GRIDCUT) {
      grid->cutoff = read_double();
    } else if (flag == COMM_SORT) {
      comm->commsortflag = read_int();
    } else if (flag == COMM_STYLE) {
      comm->commpartstyle = read_int();

    } else if (flag == NPARTICLE) {
      nparticle_file = read_bigint();
    } else if (flag == NUNSPLIT) {
      nunsplit_file = read_bigint();
    } else if (flag == NSPLIT) {
      nsplit_file = read_int();
    } else if (flag == NSUB) {
      nsub_file = read_int();
    } else if (flag == NPOINT) {
      npoint_file = read_int();
    } else if (flag == NSURF) {
      nsurf_file = read_int();

    } else error->all(FLERR,"Invalid flag in header section of restart file");

    flag = read_int();
  }
}

/* ---------------------------------------------------------------------- */

void ReadRestart::box_params()
{
  // read flags and fields until flag = -1

  int flag = read_int();
  while (flag >= 0) {

    if (flag == DIMENSION) {
      domain->dimension = read_int();
    } else if (flag == AXISYMMETRIC) {
      domain->axisymmetric = read_int();
    } else if (flag == BOXLO) {
      read_int();
      read_double_vec(3,domain->boxlo);
    } else if (flag == BOXHI) {
      read_int();
      read_double_vec(3,domain->boxhi);
    } else if (flag == BFLAG) {
      read_int();
      read_int_vec(6,domain->bflag);

    } else error->all(FLERR,"Invalid flag in header section of restart file");

    flag = read_int();
  }

  domain->print_box("Read ");
  domain->set_global_box();
}

/* ---------------------------------------------------------------------- */

void ReadRestart::particle_params()
{
  int flag = read_int();
  if (flag != PARTICLE) 
    error->all(FLERR,"Invalid flag in particle section of restart file");
  read_int();
  particle->read_restart(fp);
}

/* ---------------------------------------------------------------------- */

void ReadRestart::grid_params()
{
  int flag = read_int();
  if (flag != GRID) 
    error->all(FLERR,"Invalid flag in grid section of restart file");
  read_int();
  grid->read_restart(fp);
}

/* ---------------------------------------------------------------------- */

int ReadRestart::surf_params()
{
  int flag = read_int();
  if (flag != SURF) 
    error->all(FLERR,"Invalid flag in surf section of restart file");
  int surfexist = read_int();
  if (surfexist) surf->read_restart(fp);
  return surfexist;
}

/* ---------------------------------------------------------------------- */

void ReadRestart::file_layout()
{
  int flag = read_int();
  while (flag >= 0) {
    if (flag == MULTIPROC) {
      multiproc_file = read_int();
      if (multiproc == 0 && multiproc_file)
        error->all(FLERR,"Restart file is not a multi-proc file");
      if (multiproc && multiproc_file == 0)
        error->all(FLERR,"Restart file is a multi-proc file");

    } else error->all(FLERR,"Invalid flag in layout section of restart file");

    flag = read_int();
  }
}

/* ----------------------------------------------------------------------
   create child cells that I own
   called after Grid has stored chunk of grid cells in its restart bufs
   if skipflag = 0, all cells in restart bufs are mine
   if skipflag = 1, every Pth cell in restart bufs is mine
   cells in restart bufs are unsplit or split or sub cells
------------------------------------------------------------------------- */

void ReadRestart::create_child_cells(int skipflag)
{
  int nprocs = comm->nprocs;

  Grid::ParentCell *pcells = grid->pcells;
  int nsplit,iparent,icell,isplit,index;
  cellint id,ichild;
  double lo[3],hi[3];

  // for skipflag = 0, add all child cells in Grid restart to my Grid::cells
  // for skipflag = 1, only add every Pth cell in list

#ifdef SPARTA_MAP
  std::map<cellint,int> *hash = grid->hash;
#else
  std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

  int nlocal = grid->nlocal_restart;
  cellint *ids = grid->id_restart;
  int *nsplits = grid->nsplit_restart;

  for (int i = 0; i < nlocal; i++) {
    id = ids[i];
    nsplit = nsplits[i];

    // NOTE: need more doc of this method

    // unsplit or split cell
    // add unsplit/split cells (not sub cells) to Grid::hash as create them

    if (nsplit > 0) {
      if (skipflag && (i % nprocs != me)) continue;
      iparent = grid->id_child2parent(id,ichild);
      grid->id_child_lohi(iparent,ichild,lo,hi);
      grid->add_child_cell(id,iparent,lo,hi);
      icell = grid->nlocal - 1;
      (*hash)[id] = icell;
      grid->cells[icell].nsplit = nsplit;
      if (nsplit > 1) {
        grid->nunsplitlocal--;
        grid->add_split_cell(1);
        isplit = grid->nsplitlocal - 1;
        grid->sinfo[isplit].icell = icell;
        grid->cells[icell].isplit = isplit;
        // NOTE: need to setup csubs in SplitInfo
      }

    // sub cell
    // for skipflag, add only if I also own the corresponding split cell

    } else {
      if (skipflag && hash->find(id) == hash->end()) continue;
      index = (*hash)[id];
      grid->add_sub_cell(index,1);
      icell = grid->nlocal - 1;
      grid->cells[icell].nsplit = nsplit;
      isplit = grid->cells[icell].isplit;
      // NOTE: grid->sinfo[isplit].csubs[-nsplit] = icell;
    }
  }

  // deallocate memory in Grid

  memory->destroy(grid->id_restart);
  memory->destroy(grid->nsplit_restart);
}

/* ----------------------------------------------------------------------
   store particles that I own
   called after Particle has stored chunk of particle in its restart buf
   if skipflag = 0, all particles in restart buf are mine
   if skipflag = 1, only some particles in restart buf are mine
   use cell ID hash to convert global cell ID to local index
------------------------------------------------------------------------- */

void ReadRestart::assign_particles(int skipflag)
{
  int icell,nsplit;

#ifdef SPARTA_MAP
  std::map<cellint,int> *hash = grid->hash;
#else
  std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;

  int nlocal = particle->nlocal_restart;
  Particle::OnePartRestart *pr = particle->particle_restart;
  Particle::OnePartRestart *p;

  for (int i = 0; i < nlocal; i++) {
    p = &pr[i];
    if (skipflag && hash->find(p->icell) == hash->end()) continue;
    icell = (*hash)[p->icell];
    if (p->nsplit <= 0) 
      icell = sinfo[cells[icell].isplit].csubs[-p->nsplit];
    particle->add_particle(p->id,p->ispecies,icell,p->x,p->v,p->erot,p->ivib);
  }

  // deallocate memory in Particle

  memory->sfree(particle->particle_restart);
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// low-level fread methods
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ---------------------------------------------------------------------- */

void ReadRestart::magic_string()
{
  int n = strlen(MAGIC_STRING) + 1;
  char *str = new char[n];

  int count;
  if (me == 0) count = fread(str,sizeof(char),n,fp);
  MPI_Bcast(&count,1,MPI_INT,0,world);
  if (count < n) 
    error->all(FLERR,"Invalid SPARTA restart file");
  MPI_Bcast(str,n,MPI_CHAR,0,world);
  if (strcmp(str,MAGIC_STRING) != 0) 
    error->all(FLERR,"Invalid SPARTA restart file");
  delete [] str;
}

/* ---------------------------------------------------------------------- */

void ReadRestart::endian()
{
  int endian;
  if (me == 0) fread(&endian,sizeof(int),1,fp);
  MPI_Bcast(&endian,1,MPI_INT,0,world);
  if (endian == ENDIAN) return;
  if (endian == ENDIANSWAP)
    error->all(FLERR,"Restart file byte ordering is swapped");
  else error->all(FLERR,"Restart file byte ordering is not recognized");
}

/* ---------------------------------------------------------------------- */

int ReadRestart::version_numeric()
{
  int vn;
  if (me == 0) fread(&vn,sizeof(int),1,fp);
  MPI_Bcast(&vn,1,MPI_INT,0,world);
  if (vn != VERSION_NUMERIC) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   read an int from restart file and bcast it
------------------------------------------------------------------------- */

int ReadRestart::read_int()
{
  int value;
  if (me == 0) fread(&value,sizeof(int),1,fp);
  MPI_Bcast(&value,1,MPI_INT,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a bigint from restart file and bcast it
------------------------------------------------------------------------- */

bigint ReadRestart::read_bigint()
{
  bigint value;
  if (me == 0) fread(&value,sizeof(bigint),1,fp);
  MPI_Bcast(&value,1,MPI_SPARTA_BIGINT,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a double from restart file and bcast it
------------------------------------------------------------------------- */

double ReadRestart::read_double()
{
  double value;
  if (me == 0) fread(&value,sizeof(double),1,fp);
  MPI_Bcast(&value,1,MPI_DOUBLE,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a char string (including NULL) and bcast it
   str is allocated here, ptr is returned, caller must deallocate
------------------------------------------------------------------------- */

char *ReadRestart::read_string()
{
  int n;
  if (me == 0) fread(&n,sizeof(int),1,fp);
  MPI_Bcast(&n,1,MPI_INT,0,world);
  char *value = new char[n];
  if (me == 0) fread(value,sizeof(char),n,fp);
  MPI_Bcast(value,n,MPI_CHAR,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read vector of N ints from restart file and bcast them
------------------------------------------------------------------------- */

void ReadRestart::read_int_vec(int n, int *vec)
{
  if (me == 0) fread(vec,sizeof(int),n,fp);
  MPI_Bcast(vec,n,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   read vector of N doubles from restart file and bcast them
------------------------------------------------------------------------- */

void ReadRestart::read_double_vec(int n, double *vec)
{
  if (me == 0) fread(vec,sizeof(double),n,fp);
  MPI_Bcast(vec,n,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   read vector of N chars from restart file and bcast them
------------------------------------------------------------------------- */

void ReadRestart::read_char_vec(int n, char *vec)
{
  if (me == 0) fread(vec,sizeof(char),n,fp);
  MPI_Bcast(vec,n,MPI_CHAR,0,world);
}
