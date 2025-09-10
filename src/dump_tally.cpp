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
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "dump_tally.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

// customize by adding keyword

enum{INT,DOUBLE,BIGINT,UINT,BIGUINT,STRING};        // same as Dump

#define INVOKED_PER_TALLY 64
#define CHUNK 8

/* ---------------------------------------------------------------------- */

DumpTally::DumpTally(SPARTA *sparta, int narg, char **arg) :
  Dump(sparta, narg, arg)
{
  if (narg == 5) error->all(FLERR,"No dump tally attributes specified");

  clearstep = 1;

  // do not allow buffered output until change
  //   Dump::convert_string() to use ubuf like DumpTally::write_text() does
  // will require other dumps and computes/fixes to use datatype()
  //   like DumpTally and ComputeCollideTally do

  //buffer_allow = 1;
  //buffer_flag = 1;

  dimension = domain->dimension;

  nevery = atoi(arg[3]);

  // expand args if any have wildcard character "*"
  // ok to include trailing optional args,
  //   so long as they do not have "*" between square brackets
  // nfield may be shrunk below if extra optional args exist

  int expand = 0;
  char **earg;
  int nargnew = nfield = input->expand_args(narg-5,&arg[5],1,earg);

  if (earg != &arg[5]) expand = 1;

  // allocate field vectors

  pack_choice = new FnPtrPack[nfield];
  vtype = new int[nfield];
  field2index = new int[nfield];
  argindex = new int[nfield];

  // custom computes which the dump accesses

  ncompute = 0;
  id_compute = NULL;
  compute = NULL;

  // process attributes
  // ioptional = start of additional optional args in expanded args

  int ioptional = parse_fields(nfield,earg);

  if (ioptional < nfield)
    error->all(FLERR,"Invalid attribute in dump tally command");

  // noptional = # of optional args
  // reset nfield to subtract off optional args
  // reset ioptional to what it would be in original arg list

  int noptional = nfield - ioptional;
  nfield -= noptional;
  size_one = nfield;
  ioptional = narg - noptional;

  // setup format strings

  vformat = new char*[nfield];

  format_default = new char[4*nfield+1];
  format_default[0] = '\0';

  for (int i = 0; i < nfield; i++) {
    if (vtype[i] == DOUBLE) strcat(format_default,"%g ");
    else if (vtype[i] == INT) strcat(format_default,"%d ");
    else if (vtype[i] == BIGINT) strcat(format_default,BIGINT_FORMAT " ");
    else if (vtype[i] == UINT) strcat(format_default,"%u ");
    else if (vtype[i] == BIGUINT) strcat(format_default,BIGUINT_FORMAT " ");
    else if (vtype[i] == STRING) strcat(format_default,"%s ");
    vformat[i] = NULL;
  }

  format_column_user = new char*[size_one];
  for (int i = 0; i < size_one; i++) format_column_user[i] = NULL;

  // setup column string

  int n = 0;
  for (int iarg = 0; iarg < nfield; iarg++) n += strlen(earg[iarg]) + 2;
  columns = new char[n];
  columns[0] = '\0';
  for (int iarg = 0; iarg < nfield; iarg++) {
    strcat(columns,earg[iarg]);
    strcat(columns," ");
  }

  // if wildcard expansion occurred, free earg memory from expand_args()
  // wait to do this until after column string is created

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete [] earg[i];
    memory->sfree(earg);
  }
}

/* ---------------------------------------------------------------------- */

DumpTally::~DumpTally()
{
  delete [] pack_choice;
  delete [] vtype;
  delete [] field2index;
  delete [] argindex;

  for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
  memory->sfree(id_compute);
  delete [] compute;

  for (int i = 0; i < nfield; i++) delete [] vformat[i];
  delete [] vformat;

  delete [] columns;
}

/* ---------------------------------------------------------------------- */

void DumpTally::init_style()
{
  // setup function ptrs

  if (binary) header_choice = &DumpTally::header_binary;
  else header_choice = &DumpTally::header_item;

  if (binary) write_choice = &DumpTally::write_binary;
  else if (buffer_flag == 1) write_choice = &DumpTally::write_string;
  else write_choice = &DumpTally::write_text;

  // find current ptr for each compute

  int icompute;
  for (int i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0)
      error->all(FLERR,"Could not find dump tally compute ID");
    compute[i] = modify->compute[icompute];
  }

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpTally::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpTally::header_binary(bigint ndump)
{
  fwrite(&update->ntimestep,sizeof(bigint),1,fp);
  fwrite(&ndump,sizeof(bigint),1,fp);
  fwrite(domain->bflag,6*sizeof(int),1,fp);
  fwrite(&boxxlo,sizeof(double),1,fp);
  fwrite(&boxxhi,sizeof(double),1,fp);
  fwrite(&boxylo,sizeof(double),1,fp);
  fwrite(&boxyhi,sizeof(double),1,fp);
  fwrite(&boxzlo,sizeof(double),1,fp);
  fwrite(&boxzhi,sizeof(double),1,fp);
  fwrite(&nfield,sizeof(int),1,fp);
  if (multiproc) fwrite(&nclusterprocs,sizeof(int),1,fp);
  else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpTally::header_item(bigint ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF TALLIES\n");
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
  fprintf(fp,"%g %g\n",boxxlo,boxxhi);
  fprintf(fp,"%g %g\n",boxylo,boxyhi);
  fprintf(fp,"%g %g\n",boxzlo,boxzhi);
  fprintf(fp,"ITEM: TALLIES %s\n",columns);
}

/* ---------------------------------------------------------------------- */

int DumpTally::count()
{
  // invoke Computes for per-tally quantities
  // tally count for each compute must be equal

  int flag = 0;

  if (ncompute) {
    for (int i = 0; i < ncompute; i++)
      if (!(compute[i]->invoked_flag & INVOKED_PER_TALLY)) {
        compute[i]->compute_per_tally();
        compute[i]->invoked_flag |= INVOKED_PER_TALLY;

        surfint *dummy;
        int ntally_one = compute[i]->tallyinfo(dummy);
        if (i == 0) ntally = ntally_one;
        else if (ntally_one != ntally) flag = 1;
      }
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall)
    error->all(FLERR,"Tally count not the same for all computes in dump tally");

  // return tally count on this processor

  return ntally;
}

/* ---------------------------------------------------------------------- */

void DumpTally::pack()
{
  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);
}

/* ---------------------------------------------------------------------- */

void DumpTally::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpTally::write_binary(int n, double *mybuf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(mybuf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpTally::write_string(int n, double *mybuf)
{
  fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpTally::write_text(int n, double *mybuf)
{
  int i,j;

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < size_one; j++) {
      if (vtype[j] == DOUBLE) fprintf(fp,vformat[j],mybuf[m]);
      else if (vtype[j] == INT) fprintf(fp,vformat[j],(int) ubuf(mybuf[m]).i);
      else if (vtype[j] == BIGINT) fprintf(fp,vformat[j],(bigint) ubuf(mybuf[m]).i);
      else if (vtype[j] == UINT) fprintf(fp,vformat[j],(uint32_t) ubuf(mybuf[m]).i);
      else if (vtype[j] == BIGUINT) fprintf(fp,vformat[j],(uint64_t) ubuf(mybuf[m]).i);
      m++;
    }
    fprintf(fp,"\n");
  }
}

/* ---------------------------------------------------------------------- */

int DumpTally::parse_fields(int narg, char **arg)
{
  // customize by adding to if statement

  int i;
  for (int iarg = 0; iarg < narg; iarg++) {
    i = iarg;
    argindex[i] = -1;

    // compute value = c_ID
    // if no trailing [], then index = 0, else index = int between []

    if (strncmp(arg[iarg],"c_",2) == 0) {
      pack_choice[i] = &DumpTally::pack_compute;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump tally command");
        argindex[i] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[i] = 0;

      n = modify->find_compute(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump tally compute ID");
      if (modify->compute[n]->per_tally_flag == 0)
        error->all(FLERR,"Dump tally compute does not compute per-tally info");
      if (argindex[i] == 0 && modify->compute[n]->size_per_tally_cols != 0)
        error->all(FLERR,"Dump tally compute does not compute per-tally vector");
      if (argindex[i] > 0 && modify->compute[n]->size_per_tally_cols == 0)
        error->all(FLERR,
                   "Dump tally compute does not calculate per-tally array");
      if (argindex[i] > 0 &&
          argindex[i] > modify->compute[n]->size_per_tally_cols)
        error->all(FLERR,"Dump tally compute array is accessed out-of-range");

      // set vtype by querying compute column
      // if returns -1, error b/c it doesn't provide datatype() func

      vtype[i] = modify->compute[n]->datatype(argindex[i]);
      if (vtype[i] < 0)
        error->all(FLERR,"Dump tally compute ID does not provide datatypes");

      field2index[i] = add_compute(suffix);
      delete [] suffix;

    } else return iarg;
  }

  return narg;
}

/* ----------------------------------------------------------------------
   add Compute to list of Compute objects used by dump
   return index of where this Compute is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpTally::add_compute(char *id)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(id,id_compute[icompute]) == 0) break;
  if (icompute < ncompute) return icompute;

  id_compute = (char **)
    memory->srealloc(id_compute,(ncompute+1)*sizeof(char *),"dump:id_compute");
  delete [] compute;
  compute = new Compute*[ncompute+1];

  int n = strlen(id) + 1;
  id_compute[ncompute] = new char[n];
  strcpy(id_compute[ncompute],id);
  ncompute++;
  return ncompute-1;
}

/* ----------------------------------------------------------------------
   extraction of Compute results
------------------------------------------------------------------------- */

void DumpTally::pack_compute(int n)
{
  int index = argindex[n];
  Compute *c = compute[field2index[n]];

  if (index == 0) {
    double *vector = c->vector_tally;
    for (int i = 0; i < ntally; i++) {
      buf[n] = vector[i];
      n += size_one;
    }
  } else {
    index--;
    double **array = c->array_tally;
    for (int i = 0; i < ntally; i++) {
      buf[n] = array[i][index];
      n += size_one;
    }
  }
}
