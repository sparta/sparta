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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "dump_grid.h"
#include "update.h"
#include "domain.h"
#include "grid.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

// customize by adding keyword

enum{ID,PROC,XLO,YLO,ZLO,XHI,YHI,ZHI,XC,YC,ZC,VOL,
     COMPUTE,FIX,VARIABLE};
enum{INT,DOUBLE,BIGINT,STRING};        // same as Dump

#define INVOKED_PER_GRID 16
#define CHUNK 8

/* ---------------------------------------------------------------------- */

DumpGrid::DumpGrid(SPARTA *sparta, int narg, char **arg) :
  Dump(sparta, narg, arg)
{
  if (narg == 5) error->all(FLERR,"No dump grid attributes specified");

  clearstep = 1;
  buffer_allow = 1;
  buffer_flag = 1;

  dimension = domain->dimension;

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Dump grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

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

  // computes, fixes, variables which the dump accesses

  ncompute = 0;
  id_compute = NULL;
  compute = NULL;

  nfix = 0;
  id_fix = NULL;
  fix = NULL;

  nvariable = 0;
  id_variable = NULL;
  variable = NULL;
  vbuf = NULL;

  // process attributes
  // ioptional = start of additional optional args in expanded args

  int ioptional = parse_fields(nfield,earg);

  if (ioptional < nfield)
    error->all(FLERR,"Invalid attribute in dump grid command");

  // noptional = # of optional args
  // reset nfield to subtract off optional args
  // reset ioptional to what it would be in original arg list

  int noptional = nfield - ioptional;
  nfield -= noptional;
  size_one = nfield;
  ioptional = narg - noptional;

  // max length of per-grid variable vectors

  maxgrid = 0;

  // setup format strings

  vformat = new char*[nfield];

  format_default = new char[4*nfield+1];
  format_default[0] = '\0';

  for (int i = 0; i < nfield; i++) {
    if (vtype[i] == INT) strcat(format_default,"%d ");
    else if (vtype[i] == DOUBLE) strcat(format_default,"%g ");
    else if (vtype[i] == BIGINT) strcat(format_default,BIGINT_FORMAT " ");
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

  ncpart = 0;
  cpart = NULL;
}

/* ---------------------------------------------------------------------- */

DumpGrid::~DumpGrid()
{
  delete [] pack_choice;
  delete [] vtype;
  delete [] field2index;
  delete [] argindex;

  for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
  memory->sfree(id_compute);
  delete [] compute;

  for (int i = 0; i < nfix; i++) delete [] id_fix[i];
  memory->sfree(id_fix);
  delete [] fix;

  for (int i = 0; i < nvariable; i++) delete [] id_variable[i];
  memory->sfree(id_variable);
  delete [] variable;
  for (int i = 0; i < nvariable; i++) memory->destroy(vbuf[i]);
  delete [] vbuf;

  for (int i = 0; i < nfield; i++) delete [] vformat[i];
  delete [] vformat;

  delete [] columns;

  memory->destroy(cpart);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::init_style()
{
  // setup function ptrs

  if (binary) header_choice = &DumpGrid::header_binary;
  else header_choice = &DumpGrid::header_item;

  if (binary) write_choice = &DumpGrid::write_binary;
  else if (buffer_flag == 1) write_choice = &DumpGrid::write_string;
  else write_choice = &DumpGrid::write_text;

  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable

  int icompute;
  for (int i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0)
      error->all(FLERR,"Could not find dump grid compute ID");
    compute[i] = modify->compute[icompute];
  }

  int ifix;
  for (int i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all(FLERR,"Could not find dump grid fix ID");
    fix[i] = modify->fix[ifix];
    if (nevery % modify->fix[ifix]->per_grid_freq)
      error->all(FLERR,"Dump grid and fix not computed at compatible times");
  }

  int ivariable;
  for (int i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0)
      error->all(FLERR,"Could not find dump grid variable name");
    variable[i] = ivariable;
  }

  // create cpart index of owned grid cells with particles in grid group

  reset_grid_count();

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpGrid::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::header_binary(bigint ndump)
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

void DumpGrid::header_item(bigint ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF CELLS\n");
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
  fprintf(fp,"%g %g\n",boxxlo,boxxhi);
  fprintf(fp,"%g %g\n",boxylo,boxyhi);
  fprintf(fp,"%g %g\n",boxzlo,boxzhi);
  fprintf(fp,"ITEM: CELLS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

int DumpGrid::count()
{
  // grow variable vbuf arrays if needed

  int nglocal = grid->nlocal;
  if (nglocal > maxgrid) {
    maxgrid = grid->maxlocal;
    for (int i = 0; i < nvariable; i++) {
      memory->destroy(vbuf[i]);
      memory->create(vbuf[i],maxgrid,"dump:vbuf");
    }
  }

  // invoke Computes for per-grid quantities

  if (ncompute) {
    for (int i = 0; i < ncompute; i++)
      if (!(compute[i]->invoked_flag & INVOKED_PER_GRID)) {
        compute[i]->compute_per_grid();
        compute[i]->invoked_flag |= INVOKED_PER_GRID;
      }
  }

  // evaluate grid-style Variables for per-grid quantities

  if (nvariable)
    for (int i = 0; i < nvariable; i++)
      input->variable->compute_grid(variable[i],vbuf[i],1,0);

  // return # of grid cells with particles

  return ncpart;
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack()
{
  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::write_binary(int n, double *mybuf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(mybuf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::write_string(int n, double *mybuf)
{
  fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::write_text(int n, double *mybuf)
{
  int i,j;
  char str[32];

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < size_one; j++) {
      if (vtype[j] == INT)
        fprintf(fp,vformat[j],static_cast<int> (mybuf[m]));
      else if (vtype[j] == DOUBLE)
        fprintf(fp,vformat[j],mybuf[m]);
      else if (vtype[j] == BIGINT)
        fprintf(fp,vformat[j],static_cast<bigint> (mybuf[m]));
      else if (vtype[j] == STRING) {
        grid->id_num2str(static_cast<int> (mybuf[m]),str);
        fprintf(fp,vformat[j],str);
      }
      m++;
    }
    fprintf(fp,"\n");
  }
}

/* ---------------------------------------------------------------------- */

int DumpGrid::parse_fields(int narg, char **arg)
{
  // customize by adding to if statement

  int i;
  for (int iarg = 0; iarg < narg; iarg++) {
    i = iarg;
    argindex[i] = -1;

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[i] = &DumpGrid::pack_id;
      if (sizeof(cellint) == sizeof(smallint)) vtype[i] = INT;
      else vtype[i] = BIGINT;
    } else if (strcmp(arg[iarg],"idstr") == 0) {
      pack_choice[i] = &DumpGrid::pack_id;
      vtype[i] = STRING;
    } else if (strcmp(arg[iarg],"split") == 0) {
      pack_choice[i] = &DumpGrid::pack_split;
      vtype[i] = INT;
    } else if (strcmp(arg[iarg],"proc") == 0) {
      pack_choice[i] = &DumpGrid::pack_proc;
      vtype[i] = INT;

    } else if (strcmp(arg[iarg],"xlo") == 0) {
      pack_choice[i] = &DumpGrid::pack_xlo;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"ylo") == 0) {
      pack_choice[i] = &DumpGrid::pack_ylo;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"zlo") == 0) {
      if (dimension == 2)
        error->all(FLERR,"Invalid dump grid field for 2d simulation");
      pack_choice[i] = &DumpGrid::pack_zlo;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xhi") == 0) {
      pack_choice[i] = &DumpGrid::pack_xhi;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"yhi") == 0) {
      pack_choice[i] = &DumpGrid::pack_yhi;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"zhi") == 0) {
      if (dimension == 2)
        error->all(FLERR,"Invalid dump grid field for 2d simulation");
      pack_choice[i] = &DumpGrid::pack_zhi;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xc") == 0) {
      pack_choice[i] = &DumpGrid::pack_xc;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"yc") == 0) {
      pack_choice[i] = &DumpGrid::pack_yc;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"zc") == 0) {
      if (dimension == 2)
        error->all(FLERR,"Invalid dump grid field for 2d simulation");
      pack_choice[i] = &DumpGrid::pack_zc;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"vol") == 0) {
      pack_choice[i] = &DumpGrid::pack_vol;
      vtype[i] = DOUBLE;

    // compute value = c_ID
    // if no trailing [], then index = 0, else index = int between []

    } else if (strncmp(arg[iarg],"c_",2) == 0) {
      pack_choice[i] = &DumpGrid::pack_compute;
      vtype[i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump grid command");
        argindex[i] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[i] = 0;

      n = modify->find_compute(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump grid compute ID");
      if (modify->compute[n]->per_grid_flag == 0)
        error->all(FLERR,"Dump grid compute does not compute per-grid info");
      if (argindex[i] == 0 && modify->compute[n]->size_per_grid_cols != 0)
        error->all(FLERR,"Dump grid compute does not calculate "
                   "per-grid vector");
      if (argindex[i] > 0 && modify->compute[n]->size_per_grid_cols == 0)
        error->all(FLERR,"Dump grid compute does not calculate per-grid array");
      if (argindex[i] > 0 &&
          argindex[i] > modify->compute[n]->size_per_grid_cols)
        error->all(FLERR,"Dump grid compute array is accessed out-of-range");

      field2index[i] = add_compute(suffix);
      delete [] suffix;

    // fix value = f_ID
    // if no trailing [], then index = 0, else index = int between []
    // if index = 0 and compute stores array, expand to one value per column

    } else if (strncmp(arg[iarg],"f_",2) == 0) {
      pack_choice[i] = &DumpGrid::pack_fix;
      vtype[i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump grid command");
        argindex[i] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[i] = 0;

      n = modify->find_fix(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump grid fix ID");
      if (modify->fix[n]->per_grid_flag == 0)
        error->all(FLERR,"Dump grid fix does not compute per-grid info");
      if (argindex[i] == 0 && modify->fix[n]->size_per_grid_cols != 0)
        error->all(FLERR,"Dump grid fix does not calculate a per-grid vector");
      if (argindex[i] > 0 && modify->fix[n]->size_per_grid_cols == 0)
        error->all(FLERR,"Dump grid fix does not calculate per-grid array");
      if (argindex[i] > 0 && argindex[i] > modify->fix[n]->size_per_grid_cols)
        error->all(FLERR,"Dump grid fix array is accessed out-of-range");

      field2index[i] = add_fix(suffix);
      delete [] suffix;

    // variable value = v_name

    } else if (strncmp(arg[iarg],"v_",2) == 0) {
      pack_choice[i] = &DumpGrid::pack_variable;
      vtype[i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      argindex[i] = 0;

      n = input->variable->find(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump grid variable name");
      if (input->variable->grid_style(n) == 0)
        error->all(FLERR,"Dump grid variable is not grid-style variable");

      field2index[i] = add_variable(suffix);
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

int DumpGrid::add_compute(char *id)
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
   add Fix to list of Fix objects used by dump
   return index of where this Fix is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpGrid::add_fix(char *id)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,id_fix[ifix]) == 0) break;
  if (ifix < nfix) return ifix;

  id_fix = (char **)
    memory->srealloc(id_fix,(nfix+1)*sizeof(char *),"dump:id_fix");
  delete [] fix;
  fix = new Fix*[nfix+1];

  int n = strlen(id) + 1;
  id_fix[nfix] = new char[n];
  strcpy(id_fix[nfix],id);
  nfix++;
  return nfix-1;
}

/* ----------------------------------------------------------------------
   add Variable to list of Variables used by dump
   return index of where this Variable is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpGrid::add_variable(char *id)
{
  int ivariable;
  for (ivariable = 0; ivariable < nvariable; ivariable++)
    if (strcmp(id,id_variable[ivariable]) == 0) break;
  if (ivariable < nvariable) return ivariable;

  id_variable = (char **)
    memory->srealloc(id_variable,(nvariable+1)*sizeof(char *),
                     "dump:id_variable");
  delete [] variable;
  variable = new int[nvariable+1];
  delete [] vbuf;
  vbuf = new double*[nvariable+1];
  for (int i = 0; i <= nvariable; i++) vbuf[i] = NULL;

  int n = strlen(id) + 1;
  id_variable[nvariable] = new char[n];
  strcpy(id_variable[nvariable],id);
  nvariable++;
  return nvariable-1;
}

/* ----------------------------------------------------------------------
   create cpart array to index owned grid cells with particles in grid group
   called by init or by any operation which changes the grid during a run
     e.g. fix balance, fix adapt, fix ablate
------------------------------------------------------------------------- */

void DumpGrid::reset_grid_count()
{
  memory->destroy(cpart);
  int nglocal = grid->nlocal;
  memory->create(cpart,nglocal,"dump:cpart");

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  ncpart = 0;
  for (int i = 0; i < nglocal; i++) {
    if (!(cinfo[i].mask & groupbit)) continue;
    if (cells[i].nsplit > 1) continue;
    cpart[ncpart++] = i;
  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint DumpGrid::memory_usage()
{
  bigint bytes = Dump::memory_usage();
  bytes += memory->usage(cpart,grid->nlocal);
  return bytes;
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
------------------------------------------------------------------------- */

void DumpGrid::pack_compute(int n)
{
  int index = argindex[n];
  Compute *c = compute[field2index[n]];

  // if one of post_process flags is set,
  //   invoke post_process_grid() or invoke post_process_tally()
  // else extract from compute's vector_grid and array_grid directly
  // dump buf only stores values for grid cells with particles
  //   use cpart indices to extract needed subset

  if (c->post_process_grid_flag)
    c->post_process_grid(index,1,NULL,NULL,NULL,1);
  else if (c->post_process_isurf_grid_flag)
    c->post_process_isurf_grid();

  if (index == 0 || c->post_process_grid_flag) {
    double *vector = c->vector_grid;
    for (int i = 0; i < ncpart; i++) {
      buf[n] = vector[cpart[i]];
      n += size_one;
    }
  } else {
    index--;
    double **array = c->array_grid;
    for (int i = 0; i < ncpart; i++) {
      buf[n] = array[cpart[i]][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_fix(int n)
{
  double *vector = fix[field2index[n]]->vector_grid;
  double **array = fix[field2index[n]]->array_grid;
  int index = argindex[n];

  if (index == 0) {
    for (int i = 0; i < ncpart; i++) {
      buf[n] = vector[cpart[i]];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < ncpart; i++) {
      buf[n] = array[cpart[i]][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_variable(int n)
{
  double *vector = vbuf[field2index[n]];

  for (int i = 0; i < ncpart; i++) {
    buf[n] = vector[cpart[i]];
    n += size_one;
  }
}

/* ----------------------------------------------------------------------
   one method for every attribute dump grid can output
   the grid property is packed into buf starting at n with stride size_one
   customize a new attribute by adding a method
------------------------------------------------------------------------- */

void DumpGrid::pack_id(int n)
{
  Grid::ChildCell *cells = grid->cells;

  // NOTE: cellint (bigint) won't fit in double in some cases

  for (int i = 0; i < ncpart; i++) {
    buf[n] = cells[cpart[i]].id;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_split(int n)
{
  Grid::ChildCell *cells = grid->cells;

  for (int i = 0; i < ncpart; i++) {
    // convert to human readable format:
    //   split = 0: unsplit cell
    //   split = 1..N: split cell index + 1
    buf[n] = -cells[cpart[i]].nsplit + 1;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_proc(int n)
{
  Grid::ChildCell *cells = grid->cells;

  for (int i = 0; i < ncpart; i++) {
    buf[n] = cells[cpart[i]].proc;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_xlo(int n)
{
  Grid::ChildCell *cells = grid->cells;

  for (int i = 0; i < ncpart; i++) {
    buf[n] = cells[cpart[i]].lo[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_ylo(int n)
{
  Grid::ChildCell *cells = grid->cells;

  for (int i = 0; i < ncpart; i++) {
    buf[n] = cells[cpart[i]].lo[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_zlo(int n)
{
  Grid::ChildCell *cells = grid->cells;

  for (int i = 0; i < ncpart; i++) {
    buf[n] = cells[cpart[i]].lo[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_xhi(int n)
{
  Grid::ChildCell *cells = grid->cells;

  for (int i = 0; i < ncpart; i++) {
    buf[n] = cells[cpart[i]].hi[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_yhi(int n)
{
  Grid::ChildCell *cells = grid->cells;

  for (int i = 0; i < ncpart; i++) {
    buf[n] = cells[cpart[i]].hi[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_zhi(int n)
{
  Grid::ChildCell *cells = grid->cells;

  for (int i = 0; i < ncpart; i++) {
    buf[n] = cells[cpart[i]].hi[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_xc(int n)
{
  Grid::ChildCell *cells = grid->cells;

  for (int i = 0; i < ncpart; i++) {
    buf[n] = 0.5 * (cells[cpart[i]].lo[0] + cells[cpart[i]].hi[0]);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_yc(int n)
{
  Grid::ChildCell *cells = grid->cells;

  for (int i = 0; i < ncpart; i++) {
    buf[n] = 0.5 * (cells[cpart[i]].lo[1] + cells[cpart[i]].hi[1]);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_zc(int n)
{
  Grid::ChildCell *cells = grid->cells;

  for (int i = 0; i < ncpart; i++) {
    buf[n] = 0.5 * (cells[cpart[i]].lo[2] + cells[cpart[i]].hi[2]);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_vol(int n)
{
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int i = 0; i < ncpart; i++) {
    buf[n] = cinfo[cpart[i]].volume;
    n += size_one;
  }
}
