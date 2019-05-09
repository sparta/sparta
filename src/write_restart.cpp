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
     FNUM,NRHO,VSTREAM,TEMP_THERMAL,GRAVITY,
     SURFDIST,SURFIMPLICIT,SURFMAX,CELLMAX,SPLITMAX,
     GRIDCUT,COMM_SORT,COMM_STYLE,GRID_WEIGHT,
     DIMENSION,AXISYMMETRIC,BOXLO,BOXHI,BFLAG,
     NPARTICLE,NUNSPLIT,NSPLIT,NSUB,NPOINT,NSURF,
     SPECIES,MIXTURE,PARTICLE_CUSTOM,GRID,SURF,
     GRIDPARENT,GRIDCHILD,SURFS,PARTICLES,
     MULTIPROC,PROCSPERFILE,PERPROC};    // new fields added after PERPROC

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
   called from command() and directly from output within run loop
   file = final file name to write, except may contain a "%"
------------------------------------------------------------------------- */

void WriteRestart::write(char *file)
{
  if (update->mem_limit_grid_flag)
    update->global_mem_limit = grid->nlocal*sizeof(Grid::ChildCell);
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

  // if multiproc output:
  // open multiname file on each writing proc
  // write PROCSPERFILE into new file

  if (multiproc) {
    char *multiname = new char[strlen(file) + 16];
    char *ptr = strchr(file,'%');
    *ptr = '\0';
    sprintf(multiname,"%s%d%s",file,icluster,ptr+1);
    *ptr = '%';

    if (filewriter) {
      fpmulti = fopen(multiname,"wb");
      if (fpmulti == NULL) {
        char str[128];
        sprintf(str,"Cannot open restart file %s",multiname);
        error->one(FLERR,str);
      }
      write_int(PROCSPERFILE,nclusterprocs);
    }

    delete [] multiname;
  }

  // -----------------------------------------------
  // write the restart file(s)
  // -----------------------------------------------

  // write magic string, endian flag, numeric version
  // this is written to proc 0 file

  if (me == 0) {
    magic_string();
    endian();
    version_numeric();
  }

  // write header info
  // also box, grid, surf, particle, file layout info
  // this is written to proc 0 file

  bigint btmp = particle->nlocal;
  MPI_Allreduce(&btmp,&particle->nglobal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  if (me == 0) {
    header();
    box_info();
    grid_info();
    surf_info();
    particle_info();
    file_layout();
  }

  // write out parent grid, child grid, surfs, particles
  // this info goes into proc 0 file or multifiles

  write_parent_grid();
  write_child_grid();
  write_surfs();
  write_particles();

  // close all files

  if (me == 0) fclose(fp);
  if (multiproc && filewriter) fclose(fpmulti);
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
  // also simulation box, grid, surf, particle, file layout info

  bigint btmp = particle->nlocal;
  MPI_Allreduce(&btmp,&particle->nglobal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  if (me == 0) {
    header();
    box_info();
    grid_info();
    surf_info();
    particle_info();
    file_layout();
  }

  // communication buffer for my per-proc info = child grid cells and particles
  // max_size = largest buffer needed by any proc

  int grid_send_size = grid->size_restart();
  int particle_send_size = particle->size_restart();
  int send_size = grid_send_size + particle_send_size;

  int max_size = MAX(grid_send_size,update->global_mem_limit);
  max_size = MAX(max_size,sizeof(Particle::OnePartRestart));
  max_size += 128; // extra for size and ROUNDUP(ptr)

  int max_size_global;
  MPI_Allreduce(&max_size,&max_size_global,1,MPI_INT,MPI_MAX,world);
  max_size = max_size_global;

  char *buf;
  memory->create(buf,max_size,"write_restart:buf");
  memset(buf,0,max_size);

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

  int step_size = update->global_mem_limit/sizeof(Particle::OnePartRestart);

  // number of particles per pass
  // extra pass for grid

  int my_npasses = ceil((double)particle->nlocal/step_size)+1; 

  // output of one or more native files
  // filewriter = 1 = this proc writes to file
  // ping each proc in my cluster, receive its data, write data to file
  // else wait for ping from fileproc, send my data to fileproc
  
  int tmp,recv_size;
  MPI_Status status;
  MPI_Request request;
  
  if (filewriter) {
    int total_recv_size = 0;
    int npasses = 0;
    for (int iproc = 0; iproc < nclusterprocs; iproc++) {
      if (iproc) {
        MPI_Recv(&total_recv_size,1,MPI_INT,me+iproc,0,world,&status);
        MPI_Recv(&npasses,1,MPI_INT,me+iproc,0,world,&status);
      } else {
        total_recv_size = send_size;
        npasses = my_npasses;
      }

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
          else
            n = particle->pack_restart(&buf[n],step_size,i-1);
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
    MPI_Isend(&send_size,1,MPI_INT,fileproc,0,world,&request);
    MPI_Isend(&my_npasses,1,MPI_INT,fileproc,0,world,&request);
    for (int i = 0; i < my_npasses; i++) {
      int n = 0;
      if (i == 0)
        n = grid->pack_restart(buf);
      else
        n = particle->pack_restart(&buf[n],step_size,i-1);
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
  write_double_vec(GRAVITY,3,update->gravity);
  write_int(SURFDIST,surf->distributed);
  write_int(SURFIMPLICIT,surf->implicit);
  write_int(SURFMAX,grid->maxsurfpercell);
  write_int(CELLMAX,grid->maxcellpersurf);
  write_int(SPLITMAX,grid->maxsplitpercell);
  write_double(GRIDCUT,grid->cutoff);
  write_int(COMM_SORT,comm->commsortflag);
  write_int(COMM_STYLE,comm->commpartstyle);
  write_int(GRID_WEIGHT,grid->cellweightflag);

  write_bigint(NPARTICLE,particle->nglobal);
  write_bigint(NUNSPLIT,grid->nunsplit);
  write_int(NSPLIT,grid->nsplit);
  write_int(NSUB,grid->nsub);
  write_bigint(NSURF,surf->nsurf);

  // -1 flag signals end of header section

  int flag = -1;
  fwrite(&flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out simulation box info
------------------------------------------------------------------------- */

void WriteRestart::box_info()
{
  write_int(DIMENSION,domain->dimension);
  write_int(AXISYMMETRIC,domain->axisymmetric);
  write_double_vec(BOXLO,3,domain->boxlo);
  write_double_vec(BOXHI,3,domain->boxhi);
  write_int_vec(BFLAG,6,domain->bflag);
}

/* ----------------------------------------------------------------------
   proc 0 writes out grid info
------------------------------------------------------------------------- */

void WriteRestart::grid_info()
{
  write_int(GRID,0);
  grid->write_restart_info(fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out surface info
------------------------------------------------------------------------- */

void WriteRestart::surf_info()
{
  write_int(SURF,surf->exist);
  surf->write_restart_info(fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out particle info
------------------------------------------------------------------------- */

void WriteRestart::particle_info()
{
  write_int(SPECIES,0);
  particle->write_restart_species(fp);
  write_int(MIXTURE,0);
  particle->write_restart_mixture(fp);
  write_int(PARTICLE_CUSTOM,0);
  particle->write_restart_custom(fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out info on how file is layed out, serial or parallel
------------------------------------------------------------------------- */

void WriteRestart::file_layout()
{
  write_int(MULTIPROC,multiproc);

  // -1 flag signals end of global section

  int flag = -1;
  fwrite(&flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   write out parent cell info
   currently proc 0 writes out all parent cell data
------------------------------------------------------------------------- */

void WriteRestart::write_parent_grid()
{
  if (me == 0) write_int(GRIDPARENT,0);
  if (me == 0) grid->write_restart_parent(fp);
}

/* ----------------------------------------------------------------------
   write out child cell info
------------------------------------------------------------------------- */

void WriteRestart::write_child_grid()
{
  if (me == 0) write_int(GRIDCHILD,0);

  // communication buffer for my child cell info
  // max_size = largest buffer needed by any proc

  int send_size = grid->size_restart_child();
  int max_size;
  MPI_Allreduce(&send_size,&max_size,1,MPI_INT,MPI_MAX,world);

  char *buf;
  memory->create(buf,max_size,"write_restart:buf");
  memset(buf,0,max_size);

  // pack my child grid cell data into buf
  // n will be send_size

  int n = grid->pack_restart(buf);

  // NOTE: abstract next section?
  //comm(send_size,buf);

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
      
      // NOTE: need to use correct fp
      write_char_vec(PERPROC,recv_size,buf);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
    MPI_Rsend(buf,send_size,MPI_CHAR,fileproc,0,world);
  }

  // clean up

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out surface element info
------------------------------------------------------------------------- */

void WriteRestart::write_surfs()
{
  if (!surf->exist) return;

  if (me == 0) write_int(SURFS,0);

  if (!surf->distributed) {
    if (me == 0) surf->write_restart_all(fp);
    return;
  }

  if (!surf->implicit) {

    // communication buffer for my child cell info
    // max_size = largest buffer needed by any proc
    
    int send_size = surf->size_restart_distributed();
    int max_size;
    MPI_Allreduce(&send_size,&max_size,1,MPI_INT,MPI_MAX,world);
    
    char *buf;
    memory->create(buf,max_size,"write_restart:buf");
    memset(buf,0,max_size);

    // pack my child grid cell data into buf
    // n will be send_size
    
    int n = surf->pack_restart_distributed(buf);

    // NOTE: abstract next section?
    //comm(send_size,buf);

    // clean up
    
    memory->destroy(buf);
    return;
  }

  // communication buffer for my child cell info
  // max_size = largest buffer needed by any proc
    
  int send_size = surf->size_restart_implicit();
  int max_size;
  MPI_Allreduce(&send_size,&max_size,1,MPI_INT,MPI_MAX,world);
    
  char *buf;
  memory->create(buf,max_size,"write_restart:buf");
  memset(buf,0,max_size);

  // pack my child grid cell data into buf
  // n will be send_size
  
  int n = surf->pack_restart_implicit(buf);
  
  // NOTE: abstract next section?
  //comm(send_size,buf);
  
  // clean up
  
  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out particle info
------------------------------------------------------------------------- */

void WriteRestart::write_particles()
{
  if (me == 0) write_int(PARTICLES,0);

  //particle->write_restart(fp);
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

void WriteRestart::write_char_vec(int flag, int n, char *vec)
{
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&n,sizeof(int),1,fp);
  fwrite(vec,sizeof(char),n,fp);
}

/* ----------------------------------------------------------------------
   write a flag and vector of N chars into restart file
------------------------------------------------------------------------- */

void WriteRestart::write_char_vec(int flag, int total, int n, char *vec)
{
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&total,sizeof(int),1,fp);
  fwrite(vec,sizeof(char),n,fp);
}

/* ----------------------------------------------------------------------
   write a vector of N chars into restart file
------------------------------------------------------------------------- */

void WriteRestart::write_char_vec(int n, char *vec)
{
  fwrite(vec,sizeof(char),n,fp);
}
