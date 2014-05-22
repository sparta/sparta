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
     GRIDPARENT,SURFFLAG,NPOINT,NLINE,NTRI,
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
  // also defines species, mixtures, grid, surfs

  header(incompatible);

  box_params();
  domain->box_exist = 1;
  particle_params();
  grid_params();
  grid->exist = 1;
  surf_params();
  surf->exist = 1;

  // read file layout info

  file_layout();

  // close header file if in multiproc mode

  if (multiproc && me == 0) fclose(fp);

  // read per-proc info
  // NOTE: what needs to go here

  int maxbuf = 0;
  double *buf = NULL;
  int m,flag;

  // input of single native file
  // nprocs_file = # of chunks in file
  // proc 0 reads a chunk and bcasts it to other procs
  // each proc unpacks the atoms, saving ones in it's sub-domain
  // check for atom in sub-domain differs for orthogonal vs triclinic box

  int n;

  if (multiproc == 0) {

    for (int iproc = 0; iproc < nprocs_file; iproc++) {
      if (read_int() != PERPROC) 
        error->all(FLERR,"Invalid flag in peratom section of restart file");

      /*
      n = read_int();
      if (n > maxbuf) {
        maxbuf = n;
        memory->destroy(buf);
        memory->create(buf,maxbuf,"read_restart:buf");
      }
      read_double_vec(n,buf);

      m = 0;
      while (m < n) {
        x = &buf[m+1];
        if (triclinic) {
          domain->x2lamda(x,lamda);
          coord = lamda;
        } else coord = x;

        if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
            coord[1] >= sublo[1] && coord[1] < subhi[1] &&
            coord[2] >= sublo[2] && coord[2] < subhi[2]) {
          m += avec->unpack_restart(&buf[m]);
        } else m += static_cast<int> (buf[m]);
      }
      */

    }

    if (me == 0) fclose(fp);
  }

  // input of multiple native files with procs <= files
  // # of files = multiproc_file
  // each proc reads a subset of files, striding by nprocs
  // each proc keeps all atoms in all perproc chunks in its files

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
        fread(buf,sizeof(double),n,fp);

        m = 0;
        //while (m < n) m += avec->unpack_restart(&buf[m]);
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
  // each proc keeps all atoms in its perproc chunks in file

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
        fread(buf,sizeof(double),n,fp);

        if (i % nclusterprocs) {
          iproc = me + (i % nclusterprocs);
          MPI_Send(&n,1,MPI_INT,iproc,0,world);
          MPI_Recv(&tmp,0,MPI_INT,iproc,0,world,&status);
          MPI_Rsend(buf,n,MPI_DOUBLE,iproc,0,world);
        }

      } else if (i % nclusterprocs == me - fileproc) {
        MPI_Recv(&n,1,MPI_INT,fileproc,0,world,&status);
        if (n > maxbuf) {
          maxbuf = n;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }
        MPI_Irecv(buf,n,MPI_DOUBLE,fileproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,fileproc,0,world);
        MPI_Wait(&request,&status);
      }

      if (i % nclusterprocs == me - fileproc) {
        m = 0;
        //while (m < n) m += avec->unpack_restart(&buf[m]);
      }
    }

    if (filereader) fclose(fp);
    MPI_Comm_free(&clustercomm);
  }

  // clean-up memory

  delete [] file;
  memory->destroy(buf);

  // for multiproc files:
  // perform irregular comm to migrate atoms to correct procs

  if (multiproc) {

    // create a temporary fix to hold and migrate extra atom info
    // necessary b/c irregular will migrate atoms
  }

  // check that all particles were assigned to procs

  /*
  bigint natoms;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " atoms\n",natoms);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " atoms\n",natoms);
  }

  if (natoms != atom->natoms)
    error->all(FLERR,"Did not assign all restart atoms correctly");

  if (me == 0) {
    if (atom->nbonds) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " bonds\n",atom->nbonds);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " bonds\n",atom->nbonds);
    }
    if (atom->nangles) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " angles\n",
                          atom->nangles);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " angles\n",
                           atom->nangles);
    }
    if (atom->ndihedrals) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " dihedrals\n",
                          atom->ndihedrals);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " dihedrals\n",
                           atom->ndihedrals);
    }
    if (atom->nimpropers) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " impropers\n",
                          atom->nimpropers);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " impropers\n",
                           atom->nimpropers);
    }
  }
  */
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
      read_double_vec(3,update->vstream);
    } else if (flag == TEMP_THERMAL) {
      update->temp_thermal = read_double();
    } else if (flag == GRAVITY) {
      read_double_vec(3,update->gravity);
    } else if (flag == SURFMAX) {
      grid->maxsurfpercell = read_int();
    } else if (flag == GRIDCUT) {
      grid->cutoff = read_double();
    } else if (flag == COMM_SORT) {
      comm->commsortflag = read_int();
    } else if (flag == COMM_STYLE) {
      comm->commpartstyle = read_int();

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
      read_double_vec(3,domain->boxlo);
    } else if (flag == BOXHI) {
      read_double_vec(3,domain->boxhi);
    } else if (flag == BFLAG) {
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
  particle->read_restart_species(fp);
  particle->read_restart_mixture(fp);
}

/* ---------------------------------------------------------------------- */

void ReadRestart::grid_params()
{
}

/* ---------------------------------------------------------------------- */

void ReadRestart::surf_params()
{
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
    }
    flag = read_int();
  }
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
   do not bcast them, caller does that if required
------------------------------------------------------------------------- */

void ReadRestart::read_int_vec(int n, int *vec)
{
  if (me == 0) fread(vec,sizeof(int),n,fp);
  MPI_Bcast(vec,n,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   read vector of N doubles from restart file and bcast them
   do not bcast them, caller does that if required
------------------------------------------------------------------------- */

void ReadRestart::read_double_vec(int n, double *vec)
{
  if (me == 0) fread(vec,sizeof(double),n,fp);
  MPI_Bcast(vec,n,MPI_DOUBLE,0,world);
}
