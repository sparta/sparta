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
#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "dump.h"
#include "domain.h"
#include "update.h"
#include "input.h"
#include "grid.h"
#include "output.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

#define BIG 1.0e20
#define IBIG 2147483647
#define EPSILON 1.0e-6

#define ONEFIELD 32
#define DELTA 1048576

enum{INT,DOUBLE,BIGINT,STRING};    // many dump files

enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain

/* ---------------------------------------------------------------------- */

Dump::Dump(SPARTA *sparta, int, char **arg) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  n = strlen(arg[4]) + 1;
  filename = new char[n];
  strcpy(filename,arg[4]);

  first_flag = 0;
  flush_flag = 1;

  format = NULL;
  format_default = NULL;

  format_line_user = NULL;
  format_float_user = NULL;
  format_int_user = NULL;
  format_bigint_user = NULL;

  clearstep = 0;
  append_flag = 0;
  buffer_allow = 0;
  buffer_flag = 0;
  padflag = 0;

  maxbuf = 0;
  buf = NULL;
  maxsbuf = 0;
  sbuf = NULL;

  // parse filename for special syntax
  // if contains '%', write one file per proc and replace % with proc-ID
  // if contains '*', write one file per timestep and replace * with timestep
  // check file suffixes
  //   if ends in .bin = binary file
  //   else if ends in .gz = gzipped text file
  //   else ASCII text file

  fp = NULL;
  singlefile_opened = 0;
  compressed = 0;
  binary = 0;
  multifile = 0;

  multiproc = 0;
  nclusterprocs = nprocs;
  filewriter = 0;
  if (me == 0) filewriter = 1;
  fileproc = 0;
  multiname = NULL;

  char *ptr;
  if ((ptr = strchr(filename,'%'))) {
    multiproc = 1;
    nclusterprocs = 1;
    filewriter = 1;
    fileproc = me;
    MPI_Comm_split(world,me,0,&clustercomm);
    multiname = new char[strlen(filename) + 16];
    *ptr = '\0';
    sprintf(multiname,"%s%d%s",filename,me,ptr+1);
    *ptr = '%';
  }

  if (strchr(filename,'*')) multifile = 1;

  char *suffix = filename + strlen(filename) - strlen(".bin");
  if (suffix > filename && strcmp(suffix,".bin") == 0) binary = 1;
  suffix = filename + strlen(filename) - strlen(".gz");
  if (suffix > filename && strcmp(suffix,".gz") == 0) compressed = 1;
}

/* ---------------------------------------------------------------------- */

Dump::~Dump()
{
  delete [] id;
  delete [] style;
  delete [] filename;
  delete [] multiname;

  delete [] format;
  delete [] format_default;
  delete [] format_line_user;
  delete [] format_float_user;
  delete [] format_int_user;
  delete [] format_bigint_user;
  for (int i = 0; i < size_one; i++) delete [] format_column_user[i];
  delete [] format_column_user;

  memory->destroy(buf);
  memory->destroy(sbuf);

  if (multiproc) MPI_Comm_free(&clustercomm);

  if (multifile == 0 && fp != NULL) {
    if (compressed) {
      if (filewriter) pclose(fp);
    } else {
      if (filewriter) fclose(fp);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::init()
{
  // format = copy of default or user-specified line format

  delete [] format;
  char *str;
  if (format_line_user) str = format_line_user;
  else str = format_default;

  int n = strlen(str) + 1;
  format = new char[n];
  strcpy(format,str);

  // tokenize the format string and add space at end of each format element
  // if user-specified int/float format exists, use it instead
  // if user-specified column format exists, use it instead
  // lo priority = line, medium priority = int/float, hi priority = column

  char *ptr;
  for (int i = 0; i < size_one; i++) {
    if (i == 0) ptr = strtok(format," \0");
    else ptr = strtok(NULL," \0");
    if (ptr == NULL) error->all(FLERR,"Dump_modify format line is too short");
    delete [] vformat[i];

    if (format_column_user[i]) {
      vformat[i] = new char[strlen(format_column_user[i]) + 2];
      strcpy(vformat[i],format_column_user[i]);
    } else if (vtype[i] == INT && format_int_user) {
      vformat[i] = new char[strlen(format_int_user) + 2];
      strcpy(vformat[i],format_int_user);
    } else if (vtype[i] == DOUBLE && format_float_user) {
      vformat[i] = new char[strlen(format_float_user) + 2];
      strcpy(vformat[i],format_float_user);
    } else if (vtype[i] == BIGINT && format_bigint_user) {
      vformat[i] = new char[strlen(format_bigint_user) + 2];
      strcpy(vformat[i],format_bigint_user);
    } else {
      vformat[i] = new char[strlen(ptr) + 2];
      strcpy(vformat[i],ptr);
    }

    vformat[i] = strcat(vformat[i]," ");
  }

  // setup boundary string

  int m = 0;
  for (int idim = 0; idim < 3; idim++) {
    for (int iside = 0; iside < 2; iside++) {
      if (domain->bflag[idim*2+iside] == PERIODIC) boundstr[m++] = 'p';
      else if (domain->bflag[idim*2+iside] == OUTFLOW) boundstr[m++] = 'o';
      else if (domain->bflag[idim*2+iside] == REFLECT) boundstr[m++] = 'r';
      else if (domain->bflag[idim*2+iside] == AXISYM) boundstr[m++] = 'a';
      else if (domain->bflag[idim*2+iside] == SURFACE) boundstr[m++] = 's';
    }
    boundstr[m++] = ' ';
  }
  boundstr[8] = '\0';

  // style-specific initialization

  init_style();
}

/* ---------------------------------------------------------------------- */

void Dump::write()
{
  // if file per timestep, open new file

  if (multifile) openfile();

  // simulation box bounds

  boxxlo = domain->boxlo[0];
  boxxhi = domain->boxhi[0];
  boxylo = domain->boxlo[1];
  boxyhi = domain->boxhi[1];
  boxzlo = domain->boxlo[2];
  boxzhi = domain->boxhi[2];

  // nme = # of dump lines this proc will contribute to dump

  nme = count();

  // ntotal = total # of dump lines in snapshot
  // nmax = max # of dump lines on any proc

  bigint bnme = nme;
  MPI_Allreduce(&bnme,&ntotal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  int nmax;
  if (multiproc != nprocs) MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  else nmax = nme;

  // write timestep header
  // for multiproc,
  //   nheader = # of lines in this file via Allreduce on clustercomm

  bigint nheader = ntotal;
  if (multiproc)
    MPI_Allreduce(&bnme,&nheader,1,MPI_SPARTA_BIGINT,MPI_SUM,clustercomm);

  if (filewriter) write_header(nheader);

  // insure buf is sized for packing and communicating
  // use nmax to insure filewriter proc can receive info from others
  // limit nmax*size_one to int since used as arg in MPI calls

  if (nmax > maxbuf) {
    if ((bigint) nmax * size_one > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }

  // pack my data into buf

  pack();

 // if buffering, convert doubles into strings
  // insure sbuf is sized for communicating
  // cannot buffer if output is to binary file

  if (buffer_flag && !binary) {
    nsme = convert_string(nme,buf);
    int nsmin,nsmax;
    MPI_Allreduce(&nsme,&nsmin,1,MPI_INT,MPI_MIN,world);
    if (nsmin < 0) error->all(FLERR,"Too much buffered per-proc info for dump");
    if (multiproc != nprocs)
      MPI_Allreduce(&nsme,&nsmax,1,MPI_INT,MPI_MAX,world);
    else nsmax = nsme;
    if (nsmax > maxsbuf) {
      maxsbuf = nsmax;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }
  }

  // filewriter = 1 = this proc writes to file
  // ping each proc in my cluster, receive its data, write data to file
  // else wait for ping from fileproc, send my data to fileproc

  int tmp,nlines,nchars;
  MPI_Status status;
  MPI_Request request;

  // comm and output buf of doubles

  if (buffer_flag == 0 || binary) {
    if (filewriter) {
      for (int iproc = 0; iproc < nclusterprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(buf,maxbuf*size_one,MPI_DOUBLE,me+iproc,0,world,&request);
          MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
          MPI_Wait(&request,&status);
          MPI_Get_count(&status,MPI_DOUBLE,&nlines);
          nlines /= size_one;
        } else nlines = nme;

        write_data(nlines,buf);
      }
      if (flush_flag) fflush(fp);

    } else {
      MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
      MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,fileproc,0,world);
    }

  // comm and output sbuf = one big string of formatted values per proc

  } else {
    if (filewriter) {
      for (int iproc = 0; iproc < nclusterprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(sbuf,maxsbuf,MPI_CHAR,me+iproc,0,world,&request);
          MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
          MPI_Wait(&request,&status);
          MPI_Get_count(&status,MPI_CHAR,&nchars);
        } else nchars = nsme;

        write_data(nchars,(double *) sbuf);
      }
      if (flush_flag) fflush(fp);

    } else {
      MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
      MPI_Rsend(sbuf,nsme,MPI_CHAR,fileproc,0,world);
    }
  }

  // if file per timestep, close file if I am filewriter

  if (multifile) {
    if (compressed) {
      if (filewriter) pclose(fp);
    } else {
      if (filewriter) fclose(fp);
    }
  }
}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or gzipped
   some derived classes override this function
------------------------------------------------------------------------- */

void Dump::openfile()
{
  // single file, already opened, so just return

  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;

  // if one file per timestep, replace '*' with current timestep

  char *filecurrent = filename;
  if (multiproc) filecurrent = multiname;

  if (multifile) {
    char *filestar = filecurrent;
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    if (padflag == 0)
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
              filestar,update->ntimestep,ptr+1);
    else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
    }
    *ptr = '*';
  }

  // each proc with filewriter = 1 opens a file

  if (filewriter) {
    if (compressed) {
#ifdef SPARTA_GZIP
      char gzip[128];
      sprintf(gzip,"gzip -6 > %s",filecurrent);
#ifdef _WIN32
      fp = _popen(gzip,"wb");
#else
      fp = popen(gzip,"w");
#endif
#else
      error->one(FLERR,"Cannot open gzipped file");
#endif
    } else if (binary) {
      fp = fopen(filecurrent,"wb");
    } else if (append_flag) {
      fp = fopen(filecurrent,"a");
    } else {
      fp = fopen(filecurrent,"w");
    }

    if (fp == NULL) error->one(FLERR,"Cannot open dump file");
  } else fp = NULL;

  // delete string with timestep replaced

  if (multifile) delete [] filecurrent;
}

/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int Dump::convert_string(int n, double *mybuf)
{
  int i,j;
  char str[32];

  int offset = 0;
  int m = 0;
  for (i = 0; i < n; i++) {
    if (offset + size_one*ONEFIELD > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }

    for (j = 0; j < size_one; j++) {
      if (vtype[j] == INT)
        offset += sprintf(&sbuf[offset],vformat[j],static_cast<int> (mybuf[m]));
      else if (vtype[j] == DOUBLE)
        offset += sprintf(&sbuf[offset],vformat[j],mybuf[m]);
      else if (vtype[j] == BIGINT)
        offset += sprintf(&sbuf[offset],vformat[j],
                          static_cast<bigint> (mybuf[m]));
      else if (vtype[j] == STRING) {
        // NOTE: this is a kludge
        // assumes any STRING field from dump particle/grid/surf
        // is a grid cell ID
        // if not, might have to move this method into child classes
        // and tailor it for each dump style
        grid->id_num2str(static_cast<int> (mybuf[m]),str);
        offset += sprintf(&sbuf[offset],vformat[j],str);
      }
      m++;
    }
    offset += sprintf(&sbuf[offset],"\n");
  }

  return offset;
}

/* ----------------------------------------------------------------------
   process params common to all dumps here
   if unknown param, call modify_param specific to the dump
------------------------------------------------------------------------- */

void Dump::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal dump_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"append") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) append_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) append_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"buffer") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) buffer_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) buffer_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      if (buffer_flag && buffer_allow == 0)
        error->all(FLERR,"Dump_modify buffer yes not allowed for this style");
      iarg += 2;

    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      int idump;
      for (idump = 0; idump < output->ndump; idump++)
        if (strcmp(id,output->dump[idump]->id) == 0) break;
      int n;
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        delete [] output->var_dump[idump];
        n = strlen(&arg[iarg+1][2]) + 1;
        output->var_dump[idump] = new char[n];
        strcpy(output->var_dump[idump],&arg[iarg+1][2]);
        n = 0;
      } else {
        n = atoi(arg[iarg+1]);
        if (n <= 0) error->all(FLERR,"Illegal dump_modify command");
      }
      output->every_dump[idump] = n;
      iarg += 2;

    } else if (strcmp(arg[iarg],"fileper") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (!multiproc)
        error->all(FLERR,"Cannot use dump_modify fileper "
                   "without % in dump file name");
      int nper = atoi(arg[iarg+1]);
      if (nper <= 0) error->all(FLERR,"Illegal dump_modify command");

      multiproc = nprocs/nper;
      if (nprocs % nper) multiproc++;
      fileproc = me/nper * nper;
      int fileprocnext = MIN(fileproc+nper,nprocs);
      nclusterprocs = fileprocnext - fileproc;
      if (me == fileproc) filewriter = 1;
      else filewriter = 0;
      int icluster = fileproc/nper;

      MPI_Comm_free(&clustercomm);
      MPI_Comm_split(world,icluster,0,&clustercomm);

      delete [] multiname;
      multiname = new char[strlen(filename) + 16];
      char *ptr = strchr(filename,'%');
      *ptr = '\0';
      sprintf(multiname,"%s%d%s",filename,icluster,ptr+1);
      *ptr = '%';
      iarg += 2;

    } else if (strcmp(arg[iarg],"first") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) first_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) first_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"flush") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) flush_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) flush_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");

      if (strcmp(arg[iarg+1],"none") == 0) {
        delete [] format_line_user;
        delete [] format_int_user;
        delete [] format_bigint_user;
        delete [] format_float_user;
        format_line_user = NULL;
        format_int_user = NULL;
        format_bigint_user = NULL;
        format_float_user = NULL;
        iarg += 2;
        continue;
      }

      if (iarg+3 > narg) error->all(FLERR,"Illegal dump_modify command");

      if (strcmp(arg[iarg+1],"line") == 0) {
        delete [] format_line_user;
        int n = strlen(arg[iarg+2]) + 1;
        format_line_user = new char[n];
        strcpy(format_line_user,arg[iarg+2]);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"int") == 0) {
        delete [] format_int_user;
        int n = strlen(arg[iarg+2]) + 1;
        format_int_user = new char[n];
        strcpy(format_int_user,arg[iarg+2]);
        delete [] format_bigint_user;
        n = strlen(format_int_user) + 8;
        format_bigint_user = new char[n];
        // replace "d" in format_int_user with bigint format specifier
        // use of &str[1] removes leading '%' from BIGINT_FORMAT string
        char *ptr = strchr(format_int_user,'d');
        if (ptr == NULL)
          error->all(FLERR,
                     "Dump_modify int format does not contain d character");
        char str[8];
        sprintf(str,"%s",BIGINT_FORMAT);
        *ptr = '\0';
        sprintf(format_bigint_user,"%s%s%s",format_int_user,&str[1],ptr+1);
        *ptr = 'd';

      } else if (strcmp(arg[iarg+1],"float") == 0) {
        delete [] format_float_user;
        int n = strlen(arg[iarg+2]) + 1;
        format_float_user = new char[n];
        strcpy(format_float_user,arg[iarg+2]);

      } else {
        int i = input->inumeric(FLERR,arg[iarg+1]) - 1;
        if (i < 0 || i >= size_one)
          error->all(FLERR,"Illegal dump_modify command");
        if (format_column_user[i]) delete [] format_column_user[i];
        int n = strlen(arg[iarg+2]) + 1;
        format_column_user[i] = new char[n];
        strcpy(format_column_user[i],arg[iarg+2]);
      }
      iarg += 3;

    } else if (strcmp(arg[iarg],"nfile") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (!multiproc)
        error->all(FLERR,"Cannot use dump_modify nfile "
                   "without % in dump file name");
      int nfile = atoi(arg[iarg+1]);
      if (nfile <= 0) error->all(FLERR,"Illegal dump_modify command");
      nfile = MIN(nfile,nprocs);

      multiproc = nfile;
      int icluster = static_cast<int> ((bigint) me * nfile/nprocs);
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

      MPI_Comm_free(&clustercomm);
      MPI_Comm_split(world,icluster,0,&clustercomm);

      delete [] multiname;
      multiname = new char[strlen(filename) + 16];
      char *ptr = strchr(filename,'%');
      *ptr = '\0';
      sprintf(multiname,"%s%d%s",filename,icluster,ptr+1);
      *ptr = '%';
      iarg += 2;

    } else if (strcmp(arg[iarg],"pad") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      padflag = atoi(arg[iarg+1]);
      if (padflag < 0) error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else {
      int n = modify_param(narg-iarg,&arg[iarg]);
      if (n == 0) error->all(FLERR,"Illegal dump_modify command");
      iarg += n;
    }
  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint Dump::memory_usage()
{
  bigint bytes = memory->usage(buf,size_one*maxbuf);
  return bytes;
}
