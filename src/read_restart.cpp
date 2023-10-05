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
#include "input.h"
#include "balance_grid.h"
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

/* ---------------------------------------------------------------------- */

ReadRestart::ReadRestart(SPARTA *spa) : Pointers(spa) {}

/* ---------------------------------------------------------------------- */

void ReadRestart::command(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal read_restart command");

  if (domain->box_exist)
    error->all(FLERR,"Cannot read_restart after simulation box is defined");

  mem_limit_flag = update->global_mem_limit > 0 ||
       (update->mem_limit_grid_flag && !grid->nlocal);

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

  if (mem_limit_flag && !multiproc)
    error->all(FLERR,"Cannot (yet) use global mem/limit without "
               "% in restart file name");

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

  // process optional args

  int gridcutflag = 0;
  int balanceflag = 0;

  if (narg > 1) {
    if (strcmp(arg[1],"gridcut") == 0) {
      if (narg < 3) error->all(FLERR,"Illegal read_restart command");
      gridcutflag = 2;
    } else if (strcmp(arg[1],"balance") == 0) {
      if (narg < 3) error->all(FLERR,"Illegal read_restart command");
      balanceflag = 2;
    } else error->all(FLERR,"Illegal read_restart command");
  }

  // preform rebalancing

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  // read magic string, endian flag, numeric version

  magic_string();
  endian();
  int incompatible = version_numeric();

  // read header info which creates simulation box
  // also read particle params: species, mixture, custom attributes
  // also read parent grid and surfs

  header(incompatible);

  box_params();
  domain->box_exist = 1;
  particle_params();
  grid_params();
  grid->exist = 1;
  surf->exist = surf_params();

  if (surf->exist) {
    if (domain->dimension == 2) surf->compute_line_normal(0);
    if (domain->dimension == 3) surf->compute_tri_normal(0);
  }

  // read file layout info

  file_layout();

  // close header file if in multiproc mode

  if (multiproc && me == 0) fclose(fp);

  // read per-proc info, grid cells and particles

  int flag,value,tmp,procmatch_check,procmatch;
  bigint n_big;
  int n;
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

  particle->exist = 1;
  procmatch_check = 0;

  if (multiproc == 0 && nprocs_file == nprocs) {

    if (me == 0) {
      for (int iproc = 0; iproc < nprocs_file; iproc++) {
        tmp = fread(&value,sizeof(int),1,fp);
        if (value != PERPROC)
          error->one(FLERR,"Invalid flag in peratom section of restart file");

        if (iproc == 0) filepos = ftell(fp);

        tmp = fread(&n_big,sizeof(bigint),1,fp);
        if (n_big > MAXSMALLINT)
          error->one(FLERR,"Restart file read buffer too large, use global mem/limit");
        n = n_big;
        if (n > maxbuf) {
          maxbuf = n;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }

        if (iproc > 0) {
          tmp = fread(buf,sizeof(char),n,fp);
          MPI_Send(&n,1,MPI_INT,iproc,0,world);
          MPI_Recv(&tmp,0,MPI_INT,iproc,0,world,&status);
          MPI_Send(buf,n,MPI_CHAR,iproc,0,world);
        } else fseek(fp,filepos+sizeof(bigint)+n,SEEK_SET);
      }

      // rewind and read my chunk

      fseek(fp,filepos,SEEK_SET);
      tmp = fread(&n_big,sizeof(bigint),1,fp);
      if (n_big > MAXSMALLINT)
        error->one(FLERR,"Restart file read buffer too large, use global mem/limit");
      n = n_big;
      tmp = fread(buf,sizeof(char),n,fp);

      fclose(fp);

    } else {
      MPI_Recv(&n,1,MPI_INT,0,0,world,&status);
      if (n > maxbuf) {
        maxbuf = n;
        memory->destroy(buf);
        memory->create(buf,maxbuf,"read_restart:buf");
      }
      tmp = 0;
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

      n_big = read_bigint();
      if (n_big > MAXSMALLINT)
        error->one(FLERR,"Restart file read buffer too large, use global mem/limit");
      n = n_big;
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

    if (me == 0) fclose(fp);
  }

  // input of multiple native files with procs <= files
  // # of files = multiproc_file
  // each proc reads a subset of files, striding by nprocs
  // each proc keeps all cells/particles in all perproc chunks in its files
  // 2 versions of this: limited memory and unlimimited memory

  else if (nprocs <= multiproc_file && !mem_limit_flag) {

    // unlimited-memory version

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

      tmp = fread(&flag,sizeof(int),1,fp);
      if (flag != PROCSPERFILE)
        error->one(FLERR,"Invalid flag in peratom section of restart file");
      int procsperfile;
      tmp = fread(&procsperfile,sizeof(int),1,fp);

      for (int i = 0; i < procsperfile; i++) {
        tmp = fread(&flag,sizeof(int),1,fp);
        if (flag != PERPROC)
          error->one(FLERR,"Invalid flag in peratom section of restart file");

        tmp = fread(&n_big,sizeof(bigint),1,fp);
        if (n_big > MAXSMALLINT)
          error->one(FLERR,"Restart file read buffer too large, use global mem/limit");
        n = n_big;
        if (n > maxbuf) {
          maxbuf = n;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }
        tmp = fread(buf,sizeof(char),n,fp);

        n = grid->unpack_restart(buf);
        create_child_cells(0);
        n += particle->unpack_restart(&buf[n]);
        assign_particles(0);
      }

      fclose(fp);
    }

    delete [] procfile;
  }

  else if (nprocs <= multiproc_file && mem_limit_flag) {

    // limited-memory version

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

      tmp = fread(&flag,sizeof(int),1,fp);
      if (flag != PROCSPERFILE)
        error->one(FLERR,"Invalid flag in peratom section of restart file");
      int procsperfile;
      tmp = fread(&procsperfile,sizeof(int),1,fp);

      int step_size,npasses;

      for (int i = 0; i < procsperfile; i++) {
        tmp = fread(&flag,sizeof(int),1,fp);
        if (flag != PERPROC)
          error->one(FLERR,"Invalid flag in peratom section of restart file");

        tmp = fread(&n_big,sizeof(bigint),1,fp);

        int grid_nlocal;
        tmp = fread(&grid_nlocal,sizeof(int),1,fp);
        fseek(fp,-sizeof(int),SEEK_CUR);
        int grid_read_size = grid->size_restart(grid_nlocal);
        bigint particle_read_size = n_big - grid_read_size;
        int particle_nlocal;
        fseek(fp,grid_read_size,SEEK_CUR);
        tmp = fread(&particle_nlocal,sizeof(int),1,fp);
        fseek(fp,-(sizeof(int)+grid_read_size),SEEK_CUR);

        if (update->mem_limit_grid_flag)
          update->set_mem_limit_grid(grid_nlocal);

        int nbytes_particle = sizeof(Particle::OnePartRestart);
        int nbytes_custom = particle->sizeof_custom();
        int nbytes = nbytes_particle + nbytes_custom;

        int maxbuf_new = MIN(particle_read_size,update->global_mem_limit);
        maxbuf_new = MAX(maxbuf_new,grid_read_size);
        maxbuf_new = MAX(maxbuf_new,nbytes);
        maxbuf_new += 128; // extra for size and ROUNDUP(ptr)
        if (maxbuf_new > maxbuf) {
          maxbuf = maxbuf_new;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }

        // number of particles per pass

        step_size = MIN(particle_nlocal,update->global_mem_limit/nbytes);

        // extra pass for grid

        if (particle_nlocal == 0) npasses = 2;
        else npasses = ceil((double)particle_nlocal/step_size)+1;

        int nlocal_restart = 0;
        bigint total_read_part = 0;
        for (int ii = 0; ii < npasses; ii++) {
          if (ii == 0)
            n = grid_read_size;
          else {
            n = step_size*nbytes;
            if (ii == 1) n += IROUNDUP(sizeof(int)); // ROUNDUP(ptr)
            if (ii == npasses-1) n = particle_read_size - total_read_part;
            total_read_part += n;
          }
          tmp = fread(buf,sizeof(char),n,fp);

          if (ii == 0) {
            grid->unpack_restart(buf);
            create_child_cells(0);
          } else {
            particle->unpack_restart(buf,nlocal_restart,step_size,ii-1);
            assign_particles(0);
          }
        }

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
  // 2 versions of this: limited memory and unlimimited memory

  else if (mem_limit_flag) {

    // limited-memory version
    // nclusterprocs = # of procs in my cluster that read from one file
    // filereader = 1 if this proc reads file, else 0
    // fileproc = ID of proc in my cluster who reads from file
    // clustercomm = MPI communicator within my cluster of procs
    // procmatch = 1 if each proc in cluster gets one chunk

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
      tmp = fread(&flag,sizeof(int),1,fp);
      if (flag != PROCSPERFILE)
        error->one(FLERR,"Invalid flag in peratom section of restart file");
      tmp = fread(&procsperfile,sizeof(int),1,fp);
    }
    MPI_Bcast(&procsperfile,1,MPI_INT,0,clustercomm);

    procmatch_check = 1;
    if (procsperfile == nclusterprocs) procmatch = 1;
    else procmatch = 0;

    int iproc;
    MPI_Status status;
    MPI_Request request;

    int step_size,npasses;

    for (int i = 0; i < procsperfile; i++) {
      if (filereader) {
        tmp = fread(&flag,sizeof(int),1,fp);
        if (flag != PERPROC)
          error->one(FLERR,"Invalid flag in peratom section of restart file");

        tmp = fread(&n_big,sizeof(bigint),1,fp);

        int grid_nlocal;
        tmp = fread(&grid_nlocal,sizeof(int),1,fp);
        fseek(fp,-sizeof(int),SEEK_CUR);
        int grid_read_size = grid->size_restart(grid_nlocal);
        bigint particle_read_size = n_big - grid_read_size;
        int particle_nlocal;
        fseek(fp,grid_read_size,SEEK_CUR);
        tmp = fread(&particle_nlocal,sizeof(int),1,fp);
        fseek(fp,-(sizeof(int)+grid_read_size),SEEK_CUR);

        if (update->mem_limit_grid_flag)
          update->set_mem_limit_grid(grid_nlocal);

        int nbytes_particle = sizeof(Particle::OnePartRestart);
        int nbytes_custom = particle->sizeof_custom();
        int nbytes = nbytes_particle + nbytes_custom;

        int maxbuf_new = MIN(particle_read_size,update->global_mem_limit);
        maxbuf_new = MAX(maxbuf_new,grid_read_size);
        maxbuf_new = MAX(maxbuf_new,nbytes);
        maxbuf_new += 128; // extra for size and ROUNDUP(ptr)
        if (maxbuf_new > maxbuf) {
          maxbuf = maxbuf_new;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }

        // number of particles per pass

        step_size = MIN(particle_nlocal,update->global_mem_limit/nbytes);

        // extra pass for grid

        if (particle_nlocal == 0) npasses = 2;
        else npasses = ceil((double)particle_nlocal/step_size)+1;

        if (i % nclusterprocs) {
          iproc = me + (i % nclusterprocs);
          MPI_Send(&npasses,1,MPI_INT,iproc,0,world);
          MPI_Send(&step_size,1,MPI_INT,iproc,0,world);
        }

        int nlocal_restart = 0;
        bigint total_read_part = 0;
        for (int ii = 0; ii < npasses; ii++) {
          if (ii == 0)
            n = grid_read_size;
          else {
            n = step_size*nbytes;
            if (ii == 1) n += IROUNDUP(sizeof(int)); // ROUNDUP(ptr)
            if (ii == npasses-1) n = particle_read_size - total_read_part;
            total_read_part += n;
          }
          tmp = fread(buf,sizeof(char),n,fp);

          if (i % nclusterprocs) {
            iproc = me + (i % nclusterprocs);
            MPI_Send(&n,1,MPI_INT,iproc,0,world);
            MPI_Recv(&tmp,0,MPI_INT,iproc,0,world,&status);
            MPI_Rsend(buf,n,MPI_CHAR,iproc,0,world);
          } else if (i % nclusterprocs == me - fileproc) {
            if (ii == 0) {
              grid->unpack_restart(buf);
              create_child_cells(0);
            } else {
              particle->unpack_restart(buf,nlocal_restart,step_size,ii-1);
              assign_particles(0);
            }
          }
        }

      } else if (i % nclusterprocs == me - fileproc) {
        MPI_Recv(&npasses,1,MPI_INT,fileproc,0,world,&status);
        MPI_Recv(&step_size,1,MPI_INT,fileproc,0,world,&status);
        int nlocal_restart = 0;
        for (int ii = 0; ii < npasses; ii++) {
          MPI_Recv(&n,1,MPI_INT,fileproc,0,world,&status);
          if (n > maxbuf) {
            maxbuf = n;
            memory->destroy(buf);
            memory->create(buf,maxbuf,"read_restart:buf");
          }
          tmp = 0;
          MPI_Irecv(buf,n,MPI_CHAR,fileproc,0,world,&request);
          MPI_Send(&tmp,0,MPI_INT,fileproc,0,world);
          MPI_Wait(&request,&status);

          if (ii == 0) {
            grid->unpack_restart(buf);
            create_child_cells(0);
          } else {
            particle->unpack_restart(buf,nlocal_restart,step_size,ii-1);
            assign_particles(0);
          }
        }
      }
    }

    if (filereader) fclose(fp);
    MPI_Comm_free(&clustercomm);

  } else {

    // unlimited-memory version
    // nclusterprocs = # of procs in my cluster that read from one file
    // filereader = 1 if this proc reads file, else 0
    // fileproc = ID of proc in my cluster who reads from file
    // clustercomm = MPI communicator within my cluster of procs
    // procmatch = 1 if each proc in cluster gets one chunk

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
      tmp = fread(&flag,sizeof(int),1,fp);
      if (flag != PROCSPERFILE)
        error->one(FLERR,"Invalid flag in peratom section of restart file");
      tmp = fread(&procsperfile,sizeof(int),1,fp);
    }
    MPI_Bcast(&procsperfile,1,MPI_INT,0,clustercomm);

    procmatch_check = 1;
    if (procsperfile == nclusterprocs) procmatch = 1;
    else procmatch = 0;

    int tmp,iproc;
    MPI_Status status;
    MPI_Request request;

    for (int i = 0; i < procsperfile; i++) {
      if (filereader) {
        tmp = fread(&flag,sizeof(int),1,fp);
        if (flag != PERPROC)
          error->one(FLERR,"Invalid flag in peratom section of restart file");

        tmp = fread(&n_big,sizeof(bigint),1,fp);
        if (n_big > MAXSMALLINT)
          error->one(FLERR,"Restart file read buffer too large, use global mem/limit");
        n = n_big;
        if (n > maxbuf) {
          maxbuf = n;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }
        tmp = fread(buf,sizeof(char),n,fp);

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
        tmp = 0;
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

  // setup the grid

  if (grid->cellweightflag) grid->weight(-1,NULL);
  grid->set_maxlevel();
  grid->setup_owned();

  // clumped decomposition is maintained (if original file had it)
  //   if nprocs_file = current nprocs
  //   and each proc ends up being receiving one chunk
  // for case of procmatch_check = 1, must verify
  //   that each cluster of procs that read one file
  //   read matching # of per-proc chunks in that file,
  //   do this check via MPI_Allreduce()

  if (nprocs_file != nprocs) grid->clumped = 0;
  else if (procmatch_check) {
    int allcheck;
    MPI_Allreduce(&procmatch,&allcheck,1,MPI_INT,MPI_MIN,world);
    if (allcheck == 0) grid->clumped = 0;
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
    if (domain->dimension == 2) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " surf lines\n",
                          surf->nsurf);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " surf lines\n",
                           surf->nsurf);
    } else {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " surf triangles\n",
                          surf->nsurf);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " surf triangles\n",
                           surf->nsurf);
    }
  }

  // invoke surf and grid methods to complete surf & grid setup
  // compute normals of lines or triangles

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  if (surf->exist) {
    surf->setup_owned();
    grid->clear_surf_restart();
    grid->surf2grid(0);
  }

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // if read restart file on different number of procs, clumped will be off
  // if grid cutoff is set, will not be able to acquire ghosts
  // if surfs exist, call to grid->set_inout() below will fail
  // two command options allow user to overcome this
  // gridcut can be reset to negative value
  // re-balance can be triggered (with no output)

  if (surf->exist && grid->cutoff >= 0.0 && grid->clumped == 0) {
    if (gridcutflag) grid->cutoff = input->numeric(FLERR,arg[gridcutflag]);
    else if (balanceflag) {
      BalanceGrid *bg = new BalanceGrid(sparta);
      bg->command(narg-balanceflag,&arg[balanceflag],0);
      delete bg;
    }
  }

  MPI_Barrier(world);
  double time4 = MPI_Wtime();

  grid->acquire_ghosts();
  grid->find_neighbors();
  grid->check_uniform();
  comm->reset_neighbors();

  MPI_Barrier(world);
  double time5 = MPI_Wtime();

  if (surf->exist) {
    grid->set_inout();
    grid->type_check();
  }

  // DEBUG
  //grid->debug();

  MPI_Barrier(world);
  double time6 = MPI_Wtime();

  double time_total = time6-time1;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  read/surf2grid/rebalance/ghost/inout "
              "percent = %g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total);
    }
    if (logfile) {
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  read/surf2grid/rebalance/ghost/inout "
              "percent = %g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
              100.0*(time6-time5)/time_total);
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

    } else if (flag == DT) {
      update->dt = read_double();
    } else if (flag == TIME) {
      update->time = read_double();
      update->time_last_update = update->ntimestep;

    } else if (flag == FSTYLE) {
      update->fstyle = read_int();
    } else if (flag == FIELD) {
      read_int();
      read_double_vec(3,update->field);
    } else if (flag == FIELDID) {
      update->fieldID = read_string();

    } else if (flag == SURFS_IMPLICIT) {
      surf->implicit = read_int();
    } else if (flag == SURFS_DISTRIBUTED) {
      surf->distributed = read_int();
    } else if (flag == SURFGRID) {
      grid->surfgrid_algorithm = read_int();
    } else if (flag == SURFMAX) {
      grid->maxsurfpercell = read_int();
    } else if (flag == SPLITMAX) {
      grid->maxsplitpercell = read_int();
    } else if (flag == GRIDCUT) {
      grid->cutoff = read_double();
    } else if (flag == GRID_WEIGHT) {
      grid->cellweightflag = read_int();
    } else if (flag == COMM_SORT) {
      comm->commsortflag = read_int();
    } else if (flag == COMM_STYLE) {
      comm->commpartstyle = read_int();
    } else if (flag == SURFTALLY) {
      surf->tally_comm = read_int();
    } else if (flag == PARTICLE_REORDER) {
      update->reorder_period = read_int();
    } else if (flag == MEMLIMIT_GRID) {
      // ignore value if already set

      if (mem_limit_flag) read_int();
      else update->mem_limit_grid_flag = read_int();
    } else if (flag == MEMLIMIT) {
      // ignore value if already set

      if (mem_limit_flag) read_int();
      else update->global_mem_limit = read_int();
    } else if (flag == NPARTICLE) {
      nparticle_file = read_bigint();
    } else if (flag == NUNSPLIT) {
      nunsplit_file = read_bigint();
    } else if (flag == NSPLIT) {
      nsplit_file = read_int();
    } else if (flag == NSUB) {
      nsub_file = read_int();
    } else if (flag == NSURF) {
      nsurf_file = read_bigint();

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

  domain->print_box("  ");
  domain->set_global_box();
}

/* ---------------------------------------------------------------------- */

void ReadRestart::particle_params()
{
  int flag = read_int();
  if (flag != SPECIES)
    error->all(FLERR,"Invalid flag in particle section of restart file");
  read_int();
  particle->read_restart_species(fp);

  flag = read_int();
  if (flag != MIXTURE)
    error->all(FLERR,"Invalid flag in particle section of restart file");
  read_int();
  particle->read_restart_mixture(fp);

  flag = read_int();
  if (flag != PARTICLE_CUSTOM)
    error->all(FLERR,"Invalid flag in particle section of restart file");
  read_int();
  particle->read_restart_custom(fp);
}

/* ---------------------------------------------------------------------- */

void ReadRestart::grid_params()
{
  int flag = read_int();
  if (flag != GRID)
    error->all(FLERR,"Invalid flag in grid section of restart file");
  read_int();
  grid->read_restart(fp);

  flag = read_int();
  if (flag != GRID_CUSTOM)
    error->all(FLERR,"Invalid flag in grid section of restart file");
  read_int();
  grid->read_restart_custom(fp);

  // error check on too many bits for cell IDs
  // could occur if restart file was written with 64-bit IDs and
  //   read by code compiled for 32-bit IDs

  int maxlevel = grid->maxlevel;
  int nbits = grid->plevels[maxlevel-1].nbits + grid->plevels[maxlevel-1].newbits;
  if (nbits > sizeof(cellint)*8) {
    char str[128];
    sprintf(str,"Hierarchical grid induces cell IDs that exceed %d bits",
            (int) sizeof(cellint)*8);
    error->all(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

int ReadRestart::surf_params()
{
  int flag = read_int();
  if (flag != SURF)
    error->all(FLERR,"Invalid flag in surf section of restart file");
  int surfexist = read_int();
  if (!surfexist) return 0;

  if (surfexist) surf->read_restart(fp);

  flag = read_int();
  if (flag != SURF_CUSTOM)
    error->all(FLERR,"Invalid flag in surf section of restart file");
  read_int();
  surf->read_restart_custom(fp);

  return 1;
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

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  int level,nsplit,icell,isplit,index;
  cellint id,ichild;
  double lo[3],hi[3];

  // for skipflag = 0, add all child cells in Grid restart to my Grid::cells
  // for skipflag = 1, only add every Pth cell in list

  Grid::MyHash *hash = grid->hash;
  int nlocal = grid->nlocal_restart;
  cellint *ids = grid->id_restart;
  int *levels = grid->level_restart;
  int *nsplits = grid->nsplit_restart;

  for (int i = 0; i < nlocal; i++) {
    id = ids[i];
    level = levels[i];
    nsplit = nsplits[i];

    // unsplit or split cell
    // for skipflag == 1, add only if I own this cell
    // add as child cell to grid->cells
    // if split cell, also add split cell to sinfo
    // add unsplit/split cells (not sub cells) to Grid::hash as create them

    if (nsplit > 0) {
      if (skipflag && (i % nprocs != me)) continue;
      grid->id_lohi(id,level,boxlo,boxhi,lo,hi);
      grid->add_child_cell(id,level,lo,hi);
      icell = grid->nlocal - 1;
      (*hash)[id] = icell;
      grid->cells[icell].nsplit = nsplit;
      if (nsplit > 1) {
        grid->nunsplitlocal--;
        grid->add_split_cell(1);
        isplit = grid->nsplitlocal - 1;
        grid->cells[icell].isplit = isplit;
        grid->sinfo[isplit].icell = icell;
        grid->sinfo[isplit].csubs = grid->csubs_request(nsplit);
      }

    // sub cell
    // for skipflag, add only if I also own the corresponding split cell
    // add as sub cell to grid->cells
    // set nsplit for new sub cell and csubs in owning cell's sinfo

    } else {
      if (skipflag && hash->find(id) == hash->end()) continue;
      index = (*hash)[id];
      grid->add_sub_cell(index,1);
      icell = grid->nlocal - 1;
      grid->cells[icell].nsplit = nsplit;
      isplit = grid->cells[icell].isplit;
      grid->sinfo[isplit].csubs[-nsplit] = icell;
    }
  }

  // deallocate memory in Grid

  memory->destroy(grid->id_restart);
  memory->destroy(grid->level_restart);
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
  int icell;

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;
  Grid::MyHash *hash = grid->hash;

  int nlocal = particle->nlocal_restart;
  char *ptr = particle->particle_restart;
  int nbytes_particle = sizeof(Particle::OnePartRestart);
  int nbytes_custom = particle->sizeof_custom();
  int nbytes = nbytes_particle + nbytes_custom;
  int ncustom = particle->ncustom;

  Particle::OnePartRestart *p;

  for (int i = 0; i < nlocal; i++) {
    p = (Particle::OnePartRestart *) ptr;
    if (skipflag && hash->find(p->icell) == hash->end()) {
      ptr += nbytes;
      continue;
    }
    icell = (*hash)[p->icell];
    if (p->nsplit <= 0)
      icell = sinfo[cells[icell].isplit].csubs[-p->nsplit];
    particle->add_particle(p->id,p->ispecies,icell,p->x,p->v,p->erot,p->evib);
    ptr += nbytes_particle;
    if (ncustom) {
      particle->unpack_custom(ptr,particle->nlocal-1);
      ptr += nbytes_custom;
    }
  }

  // deallocate restart memory in Particle

  memory->sfree(particle->particle_restart);
  particle->nlocal_restart = 0;
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
  if (me == 0) int tmp = fread(&endian,sizeof(int),1,fp);
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
  if (me == 0) int tmp = fread(&vn,sizeof(int),1,fp);
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
  if (me == 0) int tmp = fread(&value,sizeof(int),1,fp);
  MPI_Bcast(&value,1,MPI_INT,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a bigint from restart file and bcast it
------------------------------------------------------------------------- */

bigint ReadRestart::read_bigint()
{
  bigint value;
  if (me == 0) int tmp = fread(&value,sizeof(bigint),1,fp);
  MPI_Bcast(&value,1,MPI_SPARTA_BIGINT,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a double from restart file and bcast it
------------------------------------------------------------------------- */

double ReadRestart::read_double()
{
  double value;
  if (me == 0) int tmp = fread(&value,sizeof(double),1,fp);
  MPI_Bcast(&value,1,MPI_DOUBLE,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a char string (including NULL) and bcast it
   str is allocated here, ptr is returned, caller must deallocate
------------------------------------------------------------------------- */

char *ReadRestart::read_string()
{
  int n,tmp;
  if (me == 0) tmp = fread(&n,sizeof(int),1,fp);
  MPI_Bcast(&n,1,MPI_INT,0,world);
  char *value = new char[n];
  if (me == 0) tmp = fread(value,sizeof(char),n,fp);
  MPI_Bcast(value,n,MPI_CHAR,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read vector of N ints from restart file and bcast them
------------------------------------------------------------------------- */

void ReadRestart::read_int_vec(int n, int *vec)
{
  if (me == 0) int tmp = fread(vec,sizeof(int),n,fp);
  MPI_Bcast(vec,n,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   read vector of N doubles from restart file and bcast them
------------------------------------------------------------------------- */

void ReadRestart::read_double_vec(int n, double *vec)
{
  if (me == 0) int tmp = fread(vec,sizeof(double),n,fp);
  MPI_Bcast(vec,n,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   read vector of N chars from restart file and bcast them
------------------------------------------------------------------------- */

void ReadRestart::read_char_vec(bigint n, char *vec)
{
  if (me == 0) int tmp = fread(vec,sizeof(char),n,fp);
  MPI_Bcast(vec,(int)n,MPI_CHAR,0,world);
}
