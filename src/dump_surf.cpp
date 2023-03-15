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
#include "dump_surf.h"
#include "update.h"
#include "domain.h"
#include "surf.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

// customize by adding keyword

enum{INT,DOUBLE,BIGINT,STRING};        // same as Dump

enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain

#define INVOKED_PER_SURF 32
#define CHUNK 8

/* ---------------------------------------------------------------------- */

DumpSurf::DumpSurf(SPARTA *sparta, int narg, char **arg) :
  Dump(sparta, narg, arg)
{
  if (narg == 5) error->all(FLERR,"No dump surf attributes specified");

  clearstep = 1;
  buffer_allow = 1;
  buffer_flag = 1;

  dimension = domain->dimension;

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Dump surf group ID does not exist");
  groupbit = surf->bitmask[igroup];

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

  // custom props, computes, fixes, variables which the dump accesses

  ncustom = 0;
  id_custom = NULL;
  custom = NULL;

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
    error->all(FLERR,"Invalid attribute in dump surf command");

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
    if (vtype[i] == INT) strcat(format_default,"%d ");
    else if (vtype[i] == DOUBLE) strcat(format_default,"%g ");
    else if (vtype[i] == BIGINT) strcat(format_default,BIGINT_FORMAT " ");
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

  // trigger setup of list of owned surf elements belonging to surf group

  firstflag = 1;
  cglobal = clocal = NULL;
  buflocal = NULL;
}

/* ---------------------------------------------------------------------- */

DumpSurf::~DumpSurf()
{
  memory->destroy(cglobal);
  memory->destroy(clocal);
  memory->destroy(buflocal);

  delete [] pack_choice;
  delete [] vtype;
  delete [] field2index;
  delete [] argindex;

  for (int i = 0; i < ncustom; i++) delete [] id_custom[i];
  memory->sfree(id_custom);
  delete [] custom;

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
}

/* ---------------------------------------------------------------------- */

void DumpSurf::init_style()
{
  distributed = surf->distributed;
  implicit = surf->implicit;

  // setup function ptrs

  if (binary) header_choice = &DumpSurf::header_binary;
  else header_choice = &DumpSurf::header_item;

  if (binary) write_choice = &DumpSurf::write_binary;
  else if (buffer_flag == 1) write_choice = &DumpSurf::write_string;
  else write_choice = &DumpSurf::write_text;

  // check that each surf custom attribute still exists

  int icustom;
  for (int i = 0; i < ncustom; i++) {
    icustom = surf->find_custom(id_custom[i]);
    if (icustom < 0)
      error->all(FLERR,"Could not find dump surf custom attribute");
    custom[i] = icustom;
  }

  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable

  int icompute;
  for (int i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0)
      error->all(FLERR,"Could not find dump surf compute ID");
    compute[i] = modify->compute[icompute];
  }

  int ifix;
  for (int i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all(FLERR,"Could not find dump surf fix ID");
    fix[i] = modify->fix[ifix];
    if (nevery % modify->fix[ifix]->per_surf_freq)
      error->all(FLERR,"Dump surf and fix not computed at compatible times");
  }

  int ivariable;
  for (int i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0)
      error->all(FLERR,"Could not find dump surf variable name");
    variable[i] = ivariable;
  }

  // open single file, one time only

  if (multifile == 0) openfile();

  // one-time setup of lists of owned elements contributing to dump
  // NOTE: will need to recalculate, if allow addition of surf elements
  // nown = # of surf elements I own
  // nchoose = # of nown surf elements in surface group
  // cglobal[] = global indices for nchoose elements
  //             used to access lines/tris in Surf
  // clocal[] = local indices for nchoose elements
  //            used to access nown data from per-surf computes,fixes,variables

  if (!firstflag) return;
  firstflag = 0;

  Surf::Line *lines;
  Surf::Tri *tris;

  if (distributed && !implicit) lines = surf->mylines;
  else lines = surf->lines;
  if (distributed && !implicit) tris = surf->mytris;
  else tris = surf->tris;

  nown = surf->nown;
  int m;

  nchoose = 0;
  for (int i = 0; i < nown; i++) {
    if (dimension == 2) {
      if (!distributed) m = me + i*nprocs;
      else m = i;
      if (lines[m].mask & groupbit) nchoose++;
    } else {
      if (!distributed) m = me + i*nprocs;
      else m = i;
      if (tris[m].mask & groupbit) nchoose++;
    }
  }

  memory->create(cglobal,nchoose,"dump/surf:cglobal");
  memory->create(clocal,nchoose,"dump/surf:clocal");
  memory->create(buflocal,nown,"dump/surf:buflocal");

  nchoose = 0;
  for (int i = 0; i < nown; i++)
    if (dimension == 2) {
      if (!distributed) m = me + i*nprocs;
      else m = i;
      if (lines[m].mask & groupbit) {
        cglobal[nchoose] = m;
        clocal[nchoose++] = i;
      }
    } else {
      if (!distributed) m = me + i*nprocs;
      else m = i;
      if (tris[m].mask & groupbit) {
        cglobal[nchoose] = m;
        clocal[nchoose++] = i;
      }
    }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpSurf::header_binary(bigint ndump)
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

void DumpSurf::header_item(bigint ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF SURFS\n");
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
  fprintf(fp,"%g %g\n",boxxlo,boxxhi);
  fprintf(fp,"%g %g\n",boxylo,boxyhi);
  fprintf(fp,"%g %g\n",boxzlo,boxzhi);
  fprintf(fp,"ITEM: SURFS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

int DumpSurf::count()
{
  // invoke Computes for per-surf quantities

  if (ncompute) {
    for (int i = 0; i < ncompute; i++)
      if (!(compute[i]->invoked_flag & INVOKED_PER_SURF)) {
        compute[i]->compute_per_grid();
        compute[i]->invoked_flag |= INVOKED_PER_SURF;
      }
  }

  // evaluate surf-style Variables for per-surf quantities

  if (nvariable)
    for (int i = 0; i < nvariable; i++)
      input->variable->compute_surf(variable[i],vbuf[i],1,0);

  // return surfs I own that are in surface group

  return nchoose;
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack()
{
  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);
}

/* ---------------------------------------------------------------------- */

void DumpSurf::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpSurf::write_binary(int n, double *mybuf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(mybuf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpSurf::write_string(int n, double *mybuf)
{
  fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpSurf::write_text(int n, double *mybuf)
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

int DumpSurf::parse_fields(int narg, char **arg)
{
  // customize by adding to if statement

  int i;
  for (int iarg = 0; iarg < narg; iarg++) {
    i = iarg;
    argindex[i] = -1;

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[i] = &DumpSurf::pack_id;
      if (sizeof(surfint) == sizeof(smallint)) vtype[i] = INT;
      else vtype[i] = BIGINT;
    } else if (strcmp(arg[iarg],"type") == 0) {
      pack_choice[i] = &DumpSurf::pack_type;
      vtype[i] = INT;

    } else if (strcmp(arg[iarg],"v1x") == 0) {
      pack_choice[i] = &DumpSurf::pack_v1x;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"v1y") == 0) {
      pack_choice[i] = &DumpSurf::pack_v1y;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"v1z") == 0) {
      if (dimension == 2)
        error->all(FLERR,"Invalid dump surf field for 2d simulation");
      pack_choice[i] = &DumpSurf::pack_v1z;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"v2x") == 0) {
      pack_choice[i] = &DumpSurf::pack_v2x;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"v2y") == 0) {
      pack_choice[i] = &DumpSurf::pack_v2y;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"v2z") == 0) {
      if (dimension == 2)
        error->all(FLERR,"Invalid dump surf field for 2d simulation");
      pack_choice[i] = &DumpSurf::pack_v2z;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"v3x") == 0) {
      if (dimension == 2)
        error->all(FLERR,"Invalid dump surf field for 2d simulation");
      pack_choice[i] = &DumpSurf::pack_v3x;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"v3y") == 0) {
      if (dimension == 2)
        error->all(FLERR,"Invalid dump surf field for 2d simulation");
      pack_choice[i] = &DumpSurf::pack_v3y;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"v3z") == 0) {
      if (dimension == 2)
        error->all(FLERR,"Invalid dump surf field for 2d simulation");
      pack_choice[i] = &DumpSurf::pack_v3z;
      vtype[i] = DOUBLE;

   // custom surf vector or array
   // if no trailing [], then arg is set to 0, else arg is int between []

    } else if (strncmp(arg[iarg],"s_",2) == 0) {
      pack_choice[i] = &DumpSurf::pack_custom;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump surf command");
        argindex[i] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[i] = 0;

      n = surf->find_custom(suffix);
      if (n < 0)
        error->all(FLERR,"Could not find dump surf custom attribute");

      vtype[i] = surf->etype[n];
      if (argindex[i] == 0 && surf->esize[n] > 0)
        error->all(FLERR,
                   "Dump surf custom attribute does not store "
                   "per-surf vector");
      if (argindex[i] > 0 && surf->esize[n] == 0)
        error->all(FLERR,
                   "Dump surf custom attribute does not store "
                   "per-surf array");
      if (argindex[i] > 0 && argindex[i] > surf->esize[n])
        error->all(FLERR,
                   "Dump surf custom attribute is accessed out-of-range");

      field2index[i] = add_custom(suffix);
      delete [] suffix;

    // compute value = c_ID
    // if no trailing [], then index = 0, else index = int between []

    } else if (strncmp(arg[iarg],"c_",2) == 0) {
      pack_choice[i] = &DumpSurf::pack_compute;
      vtype[i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump surf command");
        argindex[i] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[i] = 0;

      n = modify->find_compute(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump surf compute ID");
      if (surf->implicit)
        error->all(FLERR,"Cannot use dump surf compute with implicit surfs");
      if (modify->compute[n]->per_surf_flag == 0)
        error->all(FLERR,"Dump surf compute does not compute per-surf info");
      if (argindex[i]== 0 && modify->compute[n]->size_per_surf_cols != 0)
        error->all(FLERR,"Dump surf compute does not compute per-surf vector");
      if (argindex[i] > 0 && modify->compute[n]->size_per_surf_cols == 0)
        error->all(FLERR,
                   "Dump surf compute does not calculate per-surf array");
      if (argindex[i] > 0 &&
          argindex[i] > modify->compute[n]->size_per_surf_cols)
        error->all(FLERR,"Dump surf compute array is accessed out-of-range");

      field2index[i] = add_compute(suffix);
      delete [] suffix;

    // fix value = f_ID
    // if no trailing [], then index = 0, else index = int between []
    // if index = 0 and compute stores array, expand to one value per column

    } else if (strncmp(arg[iarg],"f_",2) == 0) {
      pack_choice[i] = &DumpSurf::pack_fix;
      vtype[i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump surf command");
        argindex[i] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[i] = 0;

      n = modify->find_fix(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump surf fix ID");
      if (surf->implicit)
        error->all(FLERR,"Cannot use dump surf fix with implicit surfs");
      if (modify->fix[n]->per_surf_flag == 0)
        error->all(FLERR,"Dump surf fix does not compute per-surf info");
      if (argindex[i]== 0 && modify->fix[n]->size_per_surf_cols != 0)
        error->all(FLERR,"Dump surf fix does not compute per-surf vector");
      if (argindex[i] > 0 && modify->fix[n]->size_per_surf_cols == 0)
        error->all(FLERR,"Dump surf fix does not compute per-surf array");
      if (argindex[i] > 0 && argindex[i] > modify->fix[n]->size_per_surf_cols)
        error->all(FLERR,"Dump surf fix array is accessed out-of-range");

      field2index[i] = add_fix(suffix);
      delete [] suffix;

    // variable value = v_name

    } else if (strncmp(arg[iarg],"v_",2) == 0) {
      pack_choice[i] = &DumpSurf::pack_variable;
      vtype[i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      argindex[i] = 0;

      n = input->variable->find(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump surf variable name");
      if (input->variable->surf_style(n) == 0)
        error->all(FLERR,"Dump surf variable is not surf-style variable");

      field2index[i] = add_variable(suffix);
      delete [] suffix;

    } else return iarg;
  }

  return narg;
}

/* ----------------------------------------------------------------------
   add Custom ID to list of surf custom attribute IDs used by dump
   return index of where this custom attribute is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpSurf::add_custom(char *id)
{
  int icustom;
  for (icustom = 0; icustom < ncustom; icustom++)
    if (strcmp(id,id_custom[icustom]) == 0) break;
  if (icustom < ncustom) return icustom;

  id_custom = (char **)
    memory->srealloc(id_custom,(ncustom+1)*sizeof(char *),"dump:id_custom");
  delete [] custom;
  custom = new int[ncustom+1];

  int n = strlen(id) + 1;
  id_custom[ncustom] = new char[n];
  strcpy(id_custom[ncustom],id);
  ncustom++;
  return ncustom-1;
}

/* ----------------------------------------------------------------------
   add Compute to list of Compute objects used by dump
   return index of where this Compute is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpSurf::add_compute(char *id)
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

int DumpSurf::add_fix(char *id)
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

int DumpSurf::add_variable(char *id)
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

void DumpSurf::pack_compute(int n)
{
  int index = argindex[n];
  Compute *c = compute[field2index[n]];
  c->post_process_surf();

  if (index == 0) {
    double *vector = c->vector_surf;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = vector[clocal[i]];
      n += size_one;
    }
  } else {
    index--;
    double **array = c->array_surf;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = array[clocal[i]][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_fix(int n)
{
  double *vector = fix[field2index[n]]->vector_surf;
  double **array = fix[field2index[n]]->array_surf;
  int index = argindex[n];

  if (index == 0) {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = vector[clocal[i]];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = array[clocal[i]][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_variable(int n)
{
  double *vector = vbuf[field2index[n]];

  // NOTE: when add surf variables, check this logic

  for (int i = 0; i < nchoose; i++) {
    buf[n] = vector[clocal[i]];
    n += size_one;
  }
}

/* ----------------------------------------------------------------------
   extraction of custom surf attribute
------------------------------------------------------------------------- */

void DumpSurf::pack_custom(int n)
{
  int m;

  int index = custom[field2index[n]];

  // for now, custom data only allowed for explicit all
  // so custom data is nlocal in length, not nown
  // when enable distributed, commented out lines replace 2 previous ones

  if (surf->etype[index] == INT) {
    if (surf->esize[index] == 0) {
      int *vector = surf->eivec[surf->ewhich[index]];
      for (int i = 0; i < nchoose; i++) {
        m = me + i*nprocs;
        buf[n] = vector[m];
        //buf[n] = vector[clocal[i]];
        n += size_one;
      }
    } else {
      int icol = argindex[n]-1;
      int **array = surf->eiarray[surf->ewhich[index]];
      for (int i = 0; i < nchoose; i++) {
        m = me + i*nprocs;
        buf[n] = array[m][icol];
        //buf[n] = array[clocal[i]][icol];
        n += size_one;
      }
    }
  } else {
    if (surf->esize[index] == 0) {
      double *vector = surf->edvec[surf->ewhich[index]];
      for (int i = 0; i < nchoose; i++) {
        m = me + i*nprocs;
        buf[n] = vector[m];
        //buf[n] = vector[clocal[i]];
        n += size_one;
      }
    } else {
      int icol = argindex[n]-1;
      double **array = surf->edarray[surf->ewhich[index]];
      for (int i = 0; i < nchoose; i++) {
        m = me + i*nprocs;
        buf[n] = array[m][icol];
        //buf[n] = array[clocal[i]][icol];
        n += size_one;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   one method for every attribute dump surf can output
   the surf property is packed into buf starting at n with stride size_one
   customize a new attribute by adding a method
------------------------------------------------------------------------- */

void DumpSurf::pack_id(int n)
{
  // NOTE: surfint (bigint) won't fit in double in some cases

  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed && !implicit) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = lines[cglobal[i]].id;
      n += size_one;
    }
  } else {
    Surf::Tri *tris;
    if (distributed && !implicit) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = tris[cglobal[i]].id;
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_type(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed && !implicit) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = lines[cglobal[i]].type;
      n += size_one;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed && !implicit) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = tris[cglobal[i]].type;
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v1x(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed && !implicit) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = lines[cglobal[i]].p1[0];
      n += size_one;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed && !implicit) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = tris[cglobal[i]].p1[0];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v1y(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed && !implicit) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = lines[cglobal[i]].p1[1];
      n += size_one;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed && !implicit) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = tris[cglobal[i]].p1[1];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v1z(int n)
{
  Surf::Tri *tris;
  if (distributed && !implicit) tris = surf->mytris;
  else tris = surf->tris;
  for (int i = 0; i < nchoose; i++) {
    buf[n] = tris[cglobal[i]].p1[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v2x(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed && !implicit) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = lines[cglobal[i]].p2[0];
      n += size_one;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed && !implicit) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = tris[cglobal[i]].p2[0];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v2y(int n)
{
  if (dimension == 2) {
    Surf::Line *lines;
    if (distributed && !implicit) lines = surf->mylines;
    else lines = surf->lines;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = lines[cglobal[i]].p2[1];
      n += size_one;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris;
    if (distributed && !implicit) tris = surf->mytris;
    else tris = surf->tris;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = tris[cglobal[i]].p2[1];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v2z(int n)
{
  Surf::Tri *tris;
  if (distributed && !implicit) tris = surf->mytris;
  else tris = surf->tris;
  for (int i = 0; i < nchoose; i++) {
    buf[n] = tris[cglobal[i]].p2[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v3x(int n)
{
  Surf::Tri *tris;
  if (distributed && !implicit) tris = surf->mytris;
  else tris = surf->tris;
  for (int i = 0; i < nchoose; i++) {
    buf[n] = tris[cglobal[i]].p3[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v3y(int n)
{
  Surf::Tri *tris;
  if (distributed && !implicit) tris = surf->mytris;
  else tris = surf->tris;
  for (int i = 0; i < nchoose; i++) {
    buf[n] = tris[cglobal[i]].p3[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v3z(int n)
{
  Surf::Tri *tris;
  if (distributed && !implicit) tris = surf->mytris;
  else tris = surf->tris;
  for (int i = 0; i < nchoose; i++) {
    buf[n] = tris[cglobal[i]].p3[2];
    n += size_one;
  }
}
