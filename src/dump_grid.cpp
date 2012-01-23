/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
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

using namespace DSMC_NS;

// customize by adding keyword

enum{ID,PROC,XLO,YLO,ZLO,XHI,YHI,ZHI,
     COMPUTE,FIX,VARIABLE};
enum{INT,DOUBLE};

enum{PERIODIC,OUTFLOW,SPECULAR};            // same as Domain

#define INVOKED_PER_GRID 16

/* ---------------------------------------------------------------------- */

DumpGrid::DumpGrid(DSMC *dsmc, int narg, char **arg) :
  Dump(dsmc, narg, arg)
{
  if (narg == 4) error->all(FLERR,"No dump grid arguments specified");

  clearstep = 1;

  nevery = atoi(arg[3]);

  // size_one may be shrunk below if additional optional args exist

  size_one = nfield = narg - 4;
  pack_choice = new FnPtrPack[nfield];
  vtype = new int[nfield];

  // computes, fixes, variables which the dump accesses

  memory->create(field2index,nfield,"dump:field2index");
  memory->create(argindex,nfield,"dump:argindex");

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
  // ioptional = start of additional optional args
  // only dump image style processes optional args

  ioptional = parse_fields(narg,arg);
  if (ioptional < narg)
    error->all(FLERR,"Invalid attribute in dump grid command");
  size_one = nfield = ioptional - 4;

  // setup format strings

  vformat = new char*[size_one];

  format_default = new char[3*size_one+1];
  format_default[0] = '\0';

  for (int i = 0; i < size_one; i++) {
    if (vtype[i] == INT) strcat(format_default,"%d ");
    else if (vtype[i] == DOUBLE) strcat(format_default,"%g ");
    vformat[i] = NULL;
  }

  // setup column string

  int n = 0;
  for (int iarg = 4; iarg < narg; iarg++) n += strlen(arg[iarg]) + 2;
  columns = new char[n];
  columns[0] = '\0';
  for (int iarg = 4; iarg < narg; iarg++) {
    strcat(columns,arg[iarg]);
    strcat(columns," ");
  }
}

/* ---------------------------------------------------------------------- */

DumpGrid::~DumpGrid()
{
  delete [] pack_choice;
  delete [] vtype;
  memory->destroy(field2index);
  memory->destroy(argindex);

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

  for (int i = 0; i < size_one; i++) delete [] vformat[i];
  delete [] vformat;

  delete [] columns;
}

/* ---------------------------------------------------------------------- */

void DumpGrid::init_style()
{
  delete [] format;
  char *str;
  if (format_user) str = format_user;
  else str = format_default;

  int n = strlen(str) + 1;
  format = new char[n];
  strcpy(format,str);

  // tokenize the format string and add space at end of each format element

  char *ptr;
  for (int i = 0; i < size_one; i++) {
    if (i == 0) ptr = strtok(format," \0");
    else ptr = strtok(NULL," \0");
    delete [] vformat[i];
    vformat[i] = new char[strlen(ptr) + 2];
    strcpy(vformat[i],ptr);
    vformat[i] = strcat(vformat[i]," ");
  }

  // setup boundary string

  int m = 0;
  for (int idim = 0; idim < 3; idim++) {
    for (int iside = 0; iside < 2; iside++) {
      if (domain->bflag[idim*2+iside] == OUTFLOW) boundstr[m++] = 'o';
      else if (domain->bflag[idim*2+iside] == PERIODIC) boundstr[m++] = 'p';
      else if (domain->bflag[idim*2+iside] == SPECULAR) boundstr[m++] = 's';
    }
    boundstr[m++] = ' ';
  }
  boundstr[8] = '\0';

  // setup function ptrs

  if (binary) header_choice = &DumpGrid::header_binary;
  else header_choice = &DumpGrid::header_item;

  if (binary) write_choice = &DumpGrid::write_binary;
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
  fwrite(&size_one,sizeof(int),1,fp);
  if (multiproc) {
    int one = 1;
    fwrite(&one,sizeof(int),1,fp);
  } else fwrite(&nprocs,sizeof(int),1,fp);
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
  return grid->nlocal;
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack()
{
  nglocal = grid->nlocal;
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

void DumpGrid::write_text(int n, double *mybuf)
{
  int i,j;

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < size_one; j++) {
      if (vtype[j] == INT) fprintf(fp,vformat[j],static_cast<int> (mybuf[m]));
      else if (vtype[j] == DOUBLE) fprintf(fp,vformat[j],mybuf[m]);
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
  for (int iarg = 4; iarg < narg; iarg++) {
    i = iarg-4;

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[i] = &DumpGrid::pack_id;
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
      pack_choice[i] = &DumpGrid::pack_zlo;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xhi") == 0) {
      pack_choice[i] = &DumpGrid::pack_xhi;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"yhi") == 0) {
      pack_choice[i] = &DumpGrid::pack_yhi;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"zhi") == 0) {
      pack_choice[i] = &DumpGrid::pack_zhi;
      vtype[i] = DOUBLE;

    // compute value = c_ID
    // if no trailing [], then arg is set to 0, else arg is int between []

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
      if (argindex[i] == 0 && modify->compute[n]->size_per_grid_cols > 0)
	error->all(FLERR,
		   "Dump grid compute does not calculate "
		   "per-grid vector");
      if (argindex[i] > 0 && modify->compute[n]->size_per_grid_cols == 0)
	error->all(FLERR,
		   "Dump grid compute does not calculate "
		   "per-grid array");
      if (argindex[i] > 0 && 
	  argindex[i] > modify->compute[n]->size_per_grid_cols)
	error->all(FLERR,
		   "Dump grid compute vector is accessed out-of-range");

      field2index[i] = add_compute(suffix);
      delete [] suffix;
      
    // fix value = f_ID
    // if no trailing [], then arg is set to 0, else arg is between []

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
	error->all(FLERR,"Dump grid fix does not compute "
		   "per-grid info");
      if (argindex[i] == 0 && modify->fix[n]->size_per_grid_cols > 0)
	error->all(FLERR,"Dump grid fix does not compute "
		   "per-grid vector");
      if (argindex[i] > 0 && modify->fix[n]->size_per_grid_cols == 0)
	error->all(FLERR,"Dump grid fix does not compute "
		   "per-grid array");
      if (argindex[i] > 0 && 
	  argindex[i] > modify->fix[n]->size_per_grid_cols)
	error->all(FLERR,"Dump grid fix vector is accessed out-of-range");

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
      if (input->variable->cell_style(n) == 0)
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
   extraction of Compute, Fix, Variable results
------------------------------------------------------------------------- */

void DumpGrid::pack_compute(int n)
{
  double *vector = compute[field2index[n]]->vector_cell;
  double **array = compute[field2index[n]]->array_cell;
  int index = argindex[n];

  printf("AAA %d\n",index);

  if (index == 0) {
    for (int i = 0; i < nglocal; i++) {
      buf[n] = vector[i];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < nglocal; i++) {
      buf[n] = array[i][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_fix(int n)
{
  double *vector = fix[field2index[n]]->vector_cell;
  double **array = fix[field2index[n]]->array_cell;
  int index = argindex[n];

  if (index == 0) {
    for (int i = 0; i < nglocal; i++) {
      buf[n] = vector[i];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < nglocal; i++) {
      buf[n] = array[i][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_variable(int n)
{
  double *vector = vbuf[field2index[n]];

  for (int i = 0; i < nglocal; i++) {
    buf[n] = vector[i];
    n += size_one;
  }
}

/* ----------------------------------------------------------------------
   one method for every attribute dump particle can output
   the particle property is packed into buf starting at n with stride size_one
   customize a new attribute by adding a method
------------------------------------------------------------------------- */

void DumpGrid::pack_id(int n)
{
  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;

  for (int i = 0; i < nglocal; i++) {
    buf[n] = cells[mycells[i]].id;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_proc(int n)
{
  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;

  for (int i = 0; i < nglocal; i++) {
    buf[n] = cells[mycells[i]].proc;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_xlo(int n)
{
  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;

  for (int i = 0; i < nglocal; i++) {
    buf[n] = cells[mycells[i]].lo[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_ylo(int n)
{
  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;

  for (int i = 0; i < nglocal; i++) {
    buf[n] = cells[mycells[i]].lo[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_zlo(int n)
{
  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;

  for (int i = 0; i < nglocal; i++) {
    buf[n] = cells[mycells[i]].lo[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_xhi(int n)
{
  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;

  for (int i = 0; i < nglocal; i++) {
    buf[n] = cells[mycells[i]].hi[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_yhi(int n)
{
  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;

  for (int i = 0; i < nglocal; i++) {
    buf[n] = cells[mycells[i]].hi[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_zhi(int n)
{
  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;

  for (int i = 0; i < nglocal; i++) {
    buf[n] = cells[mycells[i]].hi[2];
    n += size_one;
  }
}
