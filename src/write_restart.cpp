/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
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
#include "stdlib.h"
#include "write_restart.h"
#include "universe.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "grid.h"
#include "surf.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

// same as read_restart.cpp

#define MAGIC_STRING "SpartA RestartT"
#define ENDIAN 0x0001
#define ENDIANSWAP 0x1000
#define VERSION_NUMERIC 0

enum{VERSION,SMALLINT,CELLINT,BIGINT,
     UNITS,NTIMESTEP,NPROCS,
     FNUM,NRHO,VSTREAM,TEMP_THERMAL,FSTYLE,FIELD,FIELDID,
     SURFS_IMPLICIT,SURFS_DISTRIBUTED,SURFGRID,SURFMAX,
     SPLITMAX,GRIDCUT,GRID_WEIGHT,COMM_SORT,COMM_STYLE,
     SURFTALLY,PARTICLE_REORDER,MEMLIMIT_GRID,MEMLIMIT,
     DIMENSION,AXISYMMETRIC,BOXLO,BOXHI,BFLAG,
     NPARTICLE,NUNSPLIT,NSPLIT,NSUB,NPOINT,NSURF,
     SPECIES,MIXTURE,
     GRID,SURF,
     PARTICLE_CUSTOM,GRID_CUSTOM,SURF_CUSTOM,
     MULTIPROC,PROCSPERFILE,PERPROC,DT,TIME};    // new fields added after TIME

enum{NOFIELD,CFIELD,PFIELD,GFIELD};             // update.cpp

/* ---------------------------------------------------------------------- */

WriteRestart::WriteRestart(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  multiproc = 0;
}

/* ----------------------------------------------------------------------
   called as write_restart command in input script
------------------------------------------------------------------------- */

void WriteRestart::command(int narg, char **arg)
{
  if (!grid->exist)
    error->all(FLERR,"Cannot write restart file before grid is defined");
  if (narg < 1) error->all(FLERR,"Illegal write_restart command");

  // if filename contains a "*", replace with current timestep

  char *ptr;
  int n = strlen(arg[0]) + 16;
  char *file = new char[n];

  if ((ptr = strchr(arg[0],'*'))) {
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s",arg[0],update->ntimestep,ptr+1);
  } else strcpy(file,arg[0]);

  // check for multiproc output

  if (strchr(arg[0],'%')) multiproc = nprocs;
  else multiproc = 0;

  int mem_limit_flag = update->global_mem_limit > 0 ||
           (update->mem_limit_grid_flag && !grid->nlocal);
  if (mem_limit_flag && !multiproc)
    error->all(FLERR,"Cannot (yet) use global mem/limit without "
               "% in restart file name");

  // setup output style and process optional args
  // also called by Output class for periodic restart files

  multiproc_options(multiproc,narg-1,&arg[1]);

  // init entire system
  // this is probably not required

  if (comm->me == 0 && screen)
    fprintf(screen,"System init for write_restart ...\n");
  sparta->init();

  // write single restart file

  write(file);
  delete [] file;
}

/* ---------------------------------------------------------------------- */

void WriteRestart::multiproc_options(int multiproc_caller,
                                     int narg, char **arg)
{
  multiproc = multiproc_caller;

  // defaults for multiproc file writing

  nclusterprocs = nprocs;
  filewriter = 0;
  if (me == 0) filewriter = 1;
  fileproc = 0;

  if (multiproc) {
    nclusterprocs = 1;
    filewriter = 1;
    fileproc = me;
    icluster = me;
  }

  // optional args

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"fileper") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal write_restart command");
      if (!multiproc)
        error->all(FLERR,"Cannot use write_restart fileper "
                   "without % in restart file name");
      int nper = atoi(arg[iarg+1]);
      if (nper <= 0) error->all(FLERR,"Illegal write_restart command");

      multiproc = nprocs/nper;
      if (nprocs % nper) multiproc++;
      fileproc = me/nper * nper;
      int fileprocnext = MIN(fileproc+nper,nprocs);
      nclusterprocs = fileprocnext - fileproc;
      if (me == fileproc) filewriter = 1;
      else filewriter = 0;
      icluster = fileproc/nper;
      iarg += 2;

    } else if (strcmp(arg[iarg],"nfile") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal write_restart command");
      if (!multiproc)
        error->all(FLERR,"Cannot use write_restart nfile "
                   "without % in restart file name");
      int nfile = atoi(arg[iarg+1]);
      if (nfile <= 0) error->all(FLERR,"Illegal write_restart command");
      nfile = MIN(nfile,nprocs);

      multiproc = nfile;
      icluster = static_cast<int> ((bigint) me * nfile/nprocs);
      fileproc = static_cast<int> ((bigint) icluster * nprocs/nfile);
      int fcluster = static_cast<int> ((bigint) fileproc * nfile/nprocs);
      if (fcluster < icluster) fileproc++;
      int fileprocnext =
        static_cast<int> ((bigint) (icluster+1) * nprocs/nfile);
      fcluster = static_cast<int> ((bigint) fileprocnext * nfile/nprocs);
      if (fcluster < icluster+1) fileprocnext++;
      nclusterprocs = fileprocnext - fileproc;
      if (me == fileproc) filewriter = 1;
      else filewriter = 0;
      iarg += 2;

    } else error->all(FLERR,"Illegal write_restart command");
  }
}

/* ----------------------------------------------------------------------
   called from command() and directly from output within run/minimize loop
   file = final file name to write, except may contain a "%"
------------------------------------------------------------------------- */

void WriteRestart::write(char *file)
{
  if (update->mem_limit_grid_flag)
    update->set_mem_limit_grid();
  if (update->global_mem_limit > 0 ||
      (update->mem_limit_grid_flag && !grid->nlocal))
    return write_less_memory(file);

  // open single restart file or base file for multiproc case

  if (me == 0) {
    char *hfile;
    if (multiproc) {
      hfile = new char[strlen(file) + 16];
      char *ptr = strchr(file,'%');
      *ptr = '\0';
      sprintf(hfile,"%s%s%s",file,"base",ptr+1);
      *ptr = '%';
    } else hfile = file;
    fp = fopen(hfile,"wb");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open restart file %s",hfile);
      error->one(FLERR,str);
    }
    if (multiproc) delete [] hfile;
  }

  // proc 0 writes magic string, endian flag, numeric version

  if (me == 0) {
    magic_string();
    endian();
    version_numeric();
  }

  // proc 0 writes header info
  // also simulation box, particle species, grid params, surf info

  bigint btmp = particle->nlocal;
  MPI_Allreduce(&btmp,&particle->nglobal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  if (me == 0) {
    header();
    box_params();
    particle_params();
    grid_params();
    surf_params();
  }

  // communication buffer for my per-proc info = child grid cells and particles
  // max_size = largest buffer needed by any proc

  bigint send_size_big = grid->size_restart();
  send_size_big += particle->size_restart_big();
  if (send_size_big > MAXSMALLINT)
    error->one(FLERR,"Restart file write buffer too large, use global mem/limit");
  int send_size = send_size_big;

  int max_size;
  MPI_Allreduce(&send_size,&max_size,1,MPI_INT,MPI_MAX,world);

  char *buf;
  memory->create(buf,max_size,"write_restart:buf");
  memset(buf,0,max_size);

  // all procs write file layout info which may include per-proc sizes

  file_layout(send_size);

  // header info is complete
  // if multiproc output:
  //   close header file, open multiname file on each writing proc,
  //   write PROCSPERFILE into new file

  if (multiproc) {
    if (me == 0) fclose(fp);

    char *multiname = new char[strlen(file) + 16];
    char *ptr = strchr(file,'%');
    *ptr = '\0';
    sprintf(multiname,"%s%d%s",file,icluster,ptr+1);
    *ptr = '%';

    if (filewriter) {
      fp = fopen(multiname,"wb");
      if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open restart file %s",multiname);
        error->one(FLERR,str);
      }
      write_int(PROCSPERFILE,nclusterprocs);
    }

    delete [] multiname;
  }

  // pack my child grid and particle data into buf

  int n = grid->pack_restart(buf);
  n += particle->pack_restart(&buf[n]);

  // output of one or more native files
  // filewriter = 1 = this proc writes to file
  // ping each proc in my cluster, receive its data, write data to file
  // else wait for ping from fileproc, send my data to fileproc

  int tmp,recv_size;
  MPI_Status status;
  MPI_Request request;

  if (filewriter) {
    for (int iproc = 0; iproc < nclusterprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,max_size,MPI_CHAR,me+iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_CHAR,&recv_size);
      } else recv_size = send_size;

      write_char_vec(PERPROC,recv_size,buf);
    }
    fclose(fp);

  } else {
    MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
    MPI_Rsend(buf,send_size,MPI_CHAR,fileproc,0,world);
  }

  // clean up

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   called from command() and directly from output within run/minimize loop
   file = final file name to write, except may contain a "%"
------------------------------------------------------------------------- */

void WriteRestart::write_less_memory(char *file)
{
  // open single restart file or base file for multiproc case

  if (me == 0) {
    char *hfile;
    if (multiproc) {
      hfile = new char[strlen(file) + 16];
      char *ptr = strchr(file,'%');
      *ptr = '\0';
      sprintf(hfile,"%s%s%s",file,"base",ptr+1);
      *ptr = '%';
    } else hfile = file;
    fp = fopen(hfile,"wb");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open restart file %s",hfile);
      error->one(FLERR,str);
    }
    if (multiproc) delete [] hfile;
  }

  // proc 0 writes magic string, endian flag, numeric version

  if (me == 0) {
    magic_string();
    endian();
    version_numeric();
  }

  // proc 0 writes header info
  // also simulation box, particle species, parent grid cells, surf info

  bigint btmp = particle->nlocal;
  MPI_Allreduce(&btmp,&particle->nglobal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  if (me == 0) {
    header();
    box_params();
    particle_params();
    grid_params();
    surf_params();
  }

  // communication buffer for my per-proc info = child grid cells and particles
  // max_size = largest buffer needed by any proc

  int grid_send_size = grid->size_restart();
  bigint particle_send_size = particle->size_restart_big();
  bigint send_size = grid_send_size + particle_send_size;

  int nbytes_particle = sizeof(Particle::OnePartRestart);
  int nbytes_custom = particle->sizeof_custom();
  int nbytes = nbytes_particle + nbytes_custom;

  int max_size = MIN(particle_send_size,update->global_mem_limit);
  max_size = MAX(max_size,grid_send_size);
  max_size = MAX(max_size,nbytes);
  max_size += 128; // extra for size and ROUNDUP(ptr)

  int max_size_global;
  MPI_Allreduce(&max_size,&max_size_global,1,MPI_INT,MPI_MAX,world);
  max_size = max_size_global;

  char *buf;
  memory->create(buf,max_size,"write_restart:buf");
  memset(buf,0,max_size);

  // all procs write file layout info which may include per-proc sizes

  int dummy = 0;
  file_layout(dummy);

  // header info is complete
  // if multiproc output:
  //   close header file, open multiname file on each writing proc,
  //   write PROCSPERFILE into new file

  if (multiproc) {
    if (me == 0) fclose(fp);

    char *multiname = new char[strlen(file) + 16];
    char *ptr = strchr(file,'%');
    *ptr = '\0';
    sprintf(multiname,"%s%d%s",file,icluster,ptr+1);
    *ptr = '%';

    if (filewriter) {
      fp = fopen(multiname,"wb");
      if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open restart file %s",multiname);
        error->one(FLERR,str);
      }
      write_int(PROCSPERFILE,nclusterprocs);
    }

    delete [] multiname;
  }

  // pack my child grid and particle data into buf

  // number of particles per pass

  int step_size = MIN(particle->nlocal,update->global_mem_limit/nbytes);

  // extra pass for grid

  int my_npasses;
  if (particle->nlocal == 0) my_npasses = 2;
  else my_npasses = ceil((double)particle->nlocal/step_size)+1;

  // output of one or more native files
  // filewriter = 1 = this proc writes to file
  // ping each proc in my cluster, receive its data, write data to file
  // else wait for ping from fileproc, send my data to fileproc

  int tmp,recv_size;
  MPI_Status status;
  MPI_Request request;

  if (filewriter) {
    bigint total_recv_size = 0;
    int npasses = 0;
    for (int iproc = 0; iproc < nclusterprocs; iproc++) {
      if (iproc) {
        MPI_Recv(&total_recv_size,1,MPI_SPARTA_BIGINT,me+iproc,0,world,&status);
        MPI_Recv(&npasses,1,MPI_INT,me+iproc,0,world,&status);
      } else {
        total_recv_size = send_size;
        npasses = my_npasses;
      }

      bigint total_write_part = 0;
      for (int i = 0; i < npasses; i++) {
        if (iproc) {
            MPI_Irecv(buf,max_size,MPI_CHAR,me+iproc,0,world,&request);
            MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
            MPI_Wait(&request,&status);
            MPI_Get_count(&status,MPI_CHAR,&recv_size);
        } else {
          int n = 0;
          if (i == 0)
            n = grid->pack_restart(buf);
          else {
            n = step_size*nbytes;
            if (i == 1) n += IROUNDUP(sizeof(int)); // ROUNDUP(ptr)
            if (i == npasses-1) n = particle_send_size - total_write_part;
            total_write_part += n;
            particle->pack_restart(buf,step_size,i-1);
          }
          recv_size = n;
        }
        if (i == 0)
          write_char_vec(PERPROC,total_recv_size,recv_size,buf);
        else
          write_char_vec(recv_size,buf);
      }
    }
    fclose(fp);

  } else {
    bigint total_write_part = 0;
    MPI_Isend(&send_size,1,MPI_SPARTA_BIGINT,fileproc,0,world,&request);
    MPI_Isend(&my_npasses,1,MPI_INT,fileproc,0,world,&request);
    for (int i = 0; i < my_npasses; i++) {
      int n = 0;
      if (i == 0)
        n = grid->pack_restart(buf);
      else {
        n = step_size*nbytes;
        if (i == 1) n += IROUNDUP(sizeof(int)); // ROUNDUP(ptr)
        if (i == my_npasses-1) n = particle_send_size - total_write_part;
        total_write_part += n;
        particle->pack_restart(buf,step_size,i-1);
      }
      MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
      MPI_Rsend(buf,n,MPI_CHAR,fileproc,0,world);
    }
  }

  // clean up

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   proc 0 writes out problem description
------------------------------------------------------------------------- */

void WriteRestart::header()
{
  write_string(VERSION,universe->version);
  write_int(SMALLINT,sizeof(smallint));
  write_int(CELLINT,sizeof(cellint));
  write_int(BIGINT,sizeof(bigint));
  write_string(UNITS,update->unit_style);
  write_bigint(NTIMESTEP,update->ntimestep);
  write_int(NPROCS,nprocs);

  write_double(FNUM,update->fnum);
  write_double(NRHO,update->nrho);
  write_double_vec(VSTREAM,3,update->vstream);
  write_double(TEMP_THERMAL,update->temp_thermal);

  write_double(DT,update->dt);
  double time = update->time +
    (update->ntimestep - update->time_last_update) * update->dt;
  write_double(TIME,time);

  write_int(FSTYLE,update->fstyle);
  if (update->fstyle == CFIELD) write_double_vec(FIELD,3,update->field);
  else if (update->fstyle == PFIELD) write_string(FIELDID,update->fieldID);
  else if (update->fstyle == GFIELD) write_string(FIELDID,update->fieldID);

  write_int(SURFS_IMPLICIT,surf->implicit);
  write_int(SURFS_DISTRIBUTED,surf->distributed);
  write_int(SURFGRID,grid->surfgrid_algorithm);
  write_int(SURFMAX,grid->maxsurfpercell);
  write_int(SPLITMAX,grid->maxsplitpercell);
  write_double(GRIDCUT,grid->cutoff);
  write_int(GRID_WEIGHT,grid->cellweightflag);
  write_int(COMM_SORT,comm->commsortflag);
  write_int(COMM_STYLE,comm->commpartstyle);
  write_int(SURFTALLY,surf->tally_comm);
  write_int(PARTICLE_REORDER,update->reorder_period);
  write_int(MEMLIMIT_GRID,update->mem_limit_grid_flag);
  write_int(MEMLIMIT,update->global_mem_limit);

  write_bigint(NPARTICLE,particle->nglobal);
  write_bigint(NUNSPLIT,grid->nunsplit);
  write_int(NSPLIT,grid->nsplit);
  write_int(NSUB,grid->nsub);
  write_bigint(NSURF,surf->nsurf);

  // -1 flag signals end of header

  int flag = -1;
  fwrite(&flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out simulation box info
------------------------------------------------------------------------- */

void WriteRestart::box_params()
{
  write_int(DIMENSION,domain->dimension);
  write_int(AXISYMMETRIC,domain->axisymmetric);
  write_double_vec(BOXLO,3,domain->boxlo);
  write_double_vec(BOXHI,3,domain->boxhi);
  write_int_vec(BFLAG,6,domain->bflag);

  // -1 flag signals end of box info

  int flag = -1;
  fwrite(&flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out species info
------------------------------------------------------------------------- */

void WriteRestart::particle_params()
{
  write_int(SPECIES,0);
  particle->write_restart_species(fp);
  write_int(MIXTURE,0);
  particle->write_restart_mixture(fp);
  write_int(PARTICLE_CUSTOM,0);
  particle->write_restart_custom(fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out parent grid info
------------------------------------------------------------------------- */

void WriteRestart::grid_params()
{
  write_int(GRID,0);
  grid->write_restart(fp);
  write_int(GRID_CUSTOM,0);
  grid->write_restart_custom(fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out surface element into
------------------------------------------------------------------------- */

void WriteRestart::surf_params()
{
  if (!surf->exist) {
    write_int(SURF,0);
    return;
  }

  write_int(SURF,1);
  surf->write_restart(fp);
  write_int(SURF_CUSTOM,0);
  surf->write_restart_custom(fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out file layout info
   all procs call this method, only proc 0 writes to file
------------------------------------------------------------------------- */

void WriteRestart::file_layout(int)
{
  if (me == 0) write_int(MULTIPROC,multiproc);

  // -1 flag signals end of file layout info

  if (me == 0) {
    int flag = -1;
    fwrite(&flag,sizeof(int),1,fp);
  }
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// low-level fwrite methods
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ---------------------------------------------------------------------- */

void WriteRestart::magic_string()
{
  int n = strlen(MAGIC_STRING) + 1;
  char *str = new char[n];
  strcpy(str,MAGIC_STRING);
  fwrite(str,sizeof(char),n,fp);
  delete [] str;
}

/* ---------------------------------------------------------------------- */

void WriteRestart::endian()
{
  int endian = ENDIAN;
  fwrite(&endian,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void WriteRestart::version_numeric()
{
  int vn = VERSION_NUMERIC;
  fwrite(&vn,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   write a flag and an int into restart file
------------------------------------------------------------------------- */

void WriteRestart::write_int(int flag, int value)
{
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&value,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   write a flag and a bigint into restart file
------------------------------------------------------------------------- */

void WriteRestart::write_bigint(int flag, bigint value)
{
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&value,sizeof(bigint),1,fp);
}

/* ----------------------------------------------------------------------
   write a flag and a double into restart file
------------------------------------------------------------------------- */

void WriteRestart::write_double(int flag, double value)
{
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&value,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
   write a flag and a char string (including NULL) into restart file
------------------------------------------------------------------------- */

void WriteRestart::write_string(int flag, char *value)
{
  int n = strlen(value) + 1;
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&n,sizeof(int),1,fp);
  fwrite(value,sizeof(char),n,fp);
}

/* ----------------------------------------------------------------------
   write a flag and vector of N ints into restart file
------------------------------------------------------------------------- */

void WriteRestart::write_int_vec(int flag, int n, int *vec)
{
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&n,sizeof(int),1,fp);
  fwrite(vec,sizeof(int),n,fp);
}

/* ----------------------------------------------------------------------
   write a flag and vector of N doubles into restart file
------------------------------------------------------------------------- */

void WriteRestart::write_double_vec(int flag, int n, double *vec)
{
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&n,sizeof(int),1,fp);
  fwrite(vec,sizeof(double),n,fp);
}

/* ----------------------------------------------------------------------
   write a flag and vector of N chars into restart file
------------------------------------------------------------------------- */

void WriteRestart::write_char_vec(int flag, bigint n, char *vec)
{
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&n,sizeof(bigint),1,fp);
  fwrite(vec,sizeof(char),n,fp);
}

/* ----------------------------------------------------------------------
   write a flag and vector of N chars into restart file
------------------------------------------------------------------------- */

void WriteRestart::write_char_vec(int flag, bigint total, int n, char *vec)
{
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&total,sizeof(bigint),1,fp);
  fwrite(vec,sizeof(char),n,fp);
}

/* ----------------------------------------------------------------------
   write a vector of N chars into restart file
------------------------------------------------------------------------- */

void WriteRestart::write_char_vec(int n, char *vec)
{
  fwrite(vec,sizeof(char),n,fp);
}
