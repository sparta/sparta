/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "dump_molecule.h"
#include "update.h"
#include "domain.h"
#include "particle.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

// customize by adding keyword

enum{ID,TYPE,X,Y,Z,XS,YS,ZS,VX,VY,VZ,
     COMPUTE,FIX,VARIABLE};
enum{LT,LE,GT,GE,EQ,NEQ};
enum{INT,DOUBLE,STRING};

enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain

#define INVOKED_PER_MOLECULE 8

/* ---------------------------------------------------------------------- */

DumpMolecule::DumpMolecule(SPARTA *sparta, int narg, char **arg) :
  Dump(sparta, narg, arg)
{
  if (narg == 4) error->all(FLERR,"No dump molecule arguments specified");

  clearstep = 1;

  nevery = atoi(arg[3]);

  // size_one may be shrunk below if additional optional args exist

  size_one = nfield = narg - 4;
  pack_choice = new FnPtrPack[nfield];
  vtype = new int[nfield];

  nthresh = 0;
  thresh_array = NULL;
  thresh_op = NULL;
  thresh_value = NULL;

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
  if (ioptional < narg && strcmp(style,"image") != 0)
    error->all(FLERR,"Invalid attribute in dump molecule command");
  size_one = nfield = ioptional - 4;

  // molecule selection arrays

  maxlocal = 0;
  choose = NULL;
  dchoose = NULL;
  clist = NULL;

  // element names

  ntypes = particle->nspecies;
  typenames = NULL;

  // setup format strings

  vformat = new char*[size_one];

  format_default = new char[3*size_one+1];
  format_default[0] = '\0';

  for (int i = 0; i < size_one; i++) {
    if (vtype[i] == INT) strcat(format_default,"%d ");
    else if (vtype[i] == DOUBLE) strcat(format_default,"%g ");
    else if (vtype[i] == STRING) strcat(format_default,"%s ");
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

DumpMolecule::~DumpMolecule()
{
  delete [] pack_choice;
  delete [] vtype;
  memory->destroy(field2index);
  memory->destroy(argindex);

  memory->destroy(thresh_array);
  memory->destroy(thresh_op);
  memory->destroy(thresh_value);

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

  memory->destroy(choose);
  memory->destroy(dchoose);
  memory->destroy(clist);

  if (typenames) {
    for (int i = 1; i <= ntypes; i++) delete [] typenames[i];
    delete [] typenames;
  }

  for (int i = 0; i < size_one; i++) delete [] vformat[i];
  delete [] vformat;

  delete [] columns;
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::init_style()
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
      else if (domain->bflag[idim*2+iside] == REFLECT) boundstr[m++] = 'r';
      else if (domain->bflag[idim*2+iside] == AXISYM) boundstr[m++] = 'a';
      else if (domain->bflag[idim*2+iside] == SURFACE) boundstr[m++] = 's';
    }
    boundstr[m++] = ' ';
  }
  boundstr[8] = '\0';

  // setup function ptrs

  if (binary) header_choice = &DumpMolecule::header_binary;
  else header_choice = &DumpMolecule::header_item;

  if (binary) write_choice = &DumpMolecule::write_binary;
  else write_choice = &DumpMolecule::write_text;

  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable

  int icompute;
  for (int i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) 
      error->all(FLERR,"Could not find dump molecule compute ID");
    compute[i] = modify->compute[icompute];
  }

  int ifix;
  for (int i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all(FLERR,"Could not find dump molecule fix ID");
    fix[i] = modify->fix[ifix];
    if (nevery % modify->fix[ifix]->per_molecule_freq)
      error->all(FLERR,
		 "Dump molecule and fix not computed at compatible times");
  }

  int ivariable;
  for (int i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0) 
      error->all(FLERR,"Could not find dump molecule variable name");
    variable[i] = ivariable;
  }

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::header_binary(bigint ndump)
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

void DumpMolecule::header_item(bigint ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
  fprintf(fp,"%g %g\n",boxxlo,boxxhi);
  fprintf(fp,"%g %g\n",boxylo,boxyhi);
  fprintf(fp,"%g %g\n",boxzlo,boxzhi);
  fprintf(fp,"ITEM: ATOMS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

int DumpMolecule::count()
{
  int i;

  // grow choose and variable vbuf arrays if needed

  int nlocal = particle->nlocal;
  if (nlocal > maxlocal) {
    maxlocal = particle->maxlocal;

    memory->destroy(choose);
    memory->destroy(dchoose);
    memory->destroy(clist);
    memory->create(choose,maxlocal,"dump:choose");
    memory->create(dchoose,maxlocal,"dump:dchoose");
    memory->create(clist,maxlocal,"dump:clist");

    for (i = 0; i < nvariable; i++) {
      memory->destroy(vbuf[i]);
      memory->create(vbuf[i],maxlocal,"dump:vbuf");
    }
  }

  // invoke Computes for per-molecule quantities

  if (ncompute) {
    for (i = 0; i < ncompute; i++)
      if (!(compute[i]->invoked_flag & INVOKED_PER_MOLECULE)) {
	compute[i]->compute_per_molecule();
	compute[i]->invoked_flag |= INVOKED_PER_MOLECULE;
      }
  }

  // evaluate molecule-style Variables for per-molecule quantities

  if (nvariable)
    for (i = 0; i < nvariable; i++)
      input->variable->compute_molecule(variable[i],vbuf[i],1,0);

  // choose all local particles for output

  for (i = 0; i < nlocal; i++) choose[i] = 1;

  // un-choose if any threshhold criterion isn't met

  if (nthresh) {
    double *ptr;
    double value;
    int nstride;

    Particle::OnePart *particles = particle->particles;
    int nlocal = particle->nlocal;
    
    for (int ithresh = 0; ithresh < nthresh; ithresh++) {

      // customize by adding to if statement

      if (thresh_array[ithresh] == ID) {
	for (i = 0; i < nlocal; i++) dchoose[i] = particles[i].id;
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == TYPE) {
	for (i = 0; i < nlocal; i++) dchoose[i] = particles[i].ispecies + 1;
	ptr = dchoose;
	nstride = 1;

      } else if (thresh_array[ithresh] == X) {
	for (i = 0; i < nlocal; i++) dchoose[i] = particles[i].x[0];
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == Y) {
	for (i = 0; i < nlocal; i++) dchoose[i] = particles[i].x[1];
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == Z) {
	for (i = 0; i < nlocal; i++) dchoose[i] = particles[i].x[2];
	ptr = dchoose;
	nstride = 1;

      } else if (thresh_array[ithresh] == XS) {
	double boxxlo = domain->boxlo[0];
	double invxprd = 1.0/domain->xprd;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = (particles[i].x[0] - boxxlo) * invxprd;
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == YS) {
	double boxylo = domain->boxlo[1];
	double invyprd = 1.0/domain->yprd;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = (particles[i].x[1] - boxylo) * invyprd;
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == ZS) {
	double boxzlo = domain->boxlo[2];
	double invzprd = 1.0/domain->zprd;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = (particles[i].x[2] - boxzlo) * invzprd;
	ptr = dchoose;
	nstride = 1;

      } else if (thresh_array[ithresh] == VX) {
	for (i = 0; i < nlocal; i++) dchoose[i] = particles[i].v[0];
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == VY) {
	for (i = 0; i < nlocal; i++) dchoose[i] = particles[i].v[1];
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == VZ) {
	for (i = 0; i < nlocal; i++) dchoose[i] = particles[i].v[2];
	ptr = dchoose;
	nstride = 1;

      } else if (thresh_array[ithresh] == COMPUTE) {
	i = nfield + ithresh;
	if (argindex[i] == 0) {
	  ptr = compute[field2index[i]]->vector_molecule;
	  nstride = 1;
	} else {
	  ptr = &compute[field2index[i]]->array_molecule[0][argindex[i]-1];
	  nstride = compute[field2index[i]]->size_per_molecule_cols;
	}

      } else if (thresh_array[ithresh] == FIX) {
	i = nfield + ithresh;
	if (argindex[i] == 0) {
	  ptr = fix[field2index[i]]->vector_molecule;
	  nstride = 1;
	} else {
	  ptr = &fix[field2index[i]]->array_molecule[0][argindex[i]-1];
	  nstride = fix[field2index[i]]->size_per_molecule_cols;
	}

      } else if (thresh_array[ithresh] == VARIABLE) {
	i = nfield + ithresh;
	ptr = vbuf[field2index[i]];
	nstride = 1;
      }

      // unselect molecules that don't meet threshhold criterion

      value = thresh_value[ithresh];

      if (thresh_op[ithresh] == LT) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr >= value) choose[i] = 0;
      } else if (thresh_op[ithresh] == LE) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr > value) choose[i] = 0;
      } else if (thresh_op[ithresh] == GT) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr <= value) choose[i] = 0;
      } else if (thresh_op[ithresh] == GE) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr < value) choose[i] = 0;
      } else if (thresh_op[ithresh] == EQ) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr != value) choose[i] = 0;
      } else if (thresh_op[ithresh] == NEQ) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr == value) choose[i] = 0;
      }
    }
  }

  // compress choose flags into clist
  // nchoose = # of selected molecules
  // clist[i] = local index of each selected molecules

  nchoose = 0;
  for (i = 0; i < nlocal; i++)
    if (choose[i]) clist[nchoose++] = i;

  return nchoose;
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::pack()
{
  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::write_binary(int n, double *mybuf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(mybuf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::write_text(int n, double *mybuf)
{
  int i,j;

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < size_one; j++) {
      if (vtype[j] == INT) fprintf(fp,vformat[j],static_cast<int> (mybuf[m]));
      else if (vtype[j] == DOUBLE) fprintf(fp,vformat[j],mybuf[m]);
      else if (vtype[j] == STRING) 
	fprintf(fp,vformat[j],typenames[(int) mybuf[m]]);
      m++;
    }
    fprintf(fp,"\n");
  }
}

/* ---------------------------------------------------------------------- */

int DumpMolecule::parse_fields(int narg, char **arg)
{
  // customize by adding to if statement

  int i;
  for (int iarg = 4; iarg < narg; iarg++) {
    i = iarg-4;

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[i] = &DumpMolecule::pack_id;
      vtype[i] = INT;
    } else if (strcmp(arg[iarg],"type") == 0) {
      pack_choice[i] = &DumpMolecule::pack_type;
      vtype[i] = INT;

    } else if (strcmp(arg[iarg],"x") == 0) {
      pack_choice[i] = &DumpMolecule::pack_x;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"y") == 0) {
      pack_choice[i] = &DumpMolecule::pack_y;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"z") == 0) {
      pack_choice[i] = &DumpMolecule::pack_z;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xs") == 0) {
      pack_choice[i] = &DumpMolecule::pack_xs;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"ys") == 0) {
      pack_choice[i] = &DumpMolecule::pack_ys;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"zs") == 0) {
      pack_choice[i] = &DumpMolecule::pack_zs;
      vtype[i] = DOUBLE;

    } else if (strcmp(arg[iarg],"vx") == 0) {
      pack_choice[i] = &DumpMolecule::pack_vx;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      pack_choice[i] = &DumpMolecule::pack_vy;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      pack_choice[i] = &DumpMolecule::pack_vz;
      vtype[i] = DOUBLE;

    // compute value = c_ID
    // if no trailing [], then arg is set to 0, else arg is int between []

    } else if (strncmp(arg[iarg],"c_",2) == 0) {
      pack_choice[i] = &DumpMolecule::pack_compute;
      vtype[i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all(FLERR,"Invalid attribute in dump molecule command");
	argindex[i] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[i] = 0;

      n = modify->find_compute(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump molecule compute ID");
      if (modify->compute[n]->per_molecule_flag == 0)
	error->all(FLERR,
		   "Dump molecule compute does not compute per-molecule info");
      if (argindex[i] == 0 && modify->compute[n]->size_per_molecule_cols > 0)
	error->all(FLERR,
		   "Dump molecule compute does not calculate "
		   "per-molecule vector");
      if (argindex[i] > 0 && modify->compute[n]->size_per_molecule_cols == 0)
	error->all(FLERR,
		   "Dump molecule compute does not calculate "
		   "per-molecule array");
      if (argindex[i] > 0 && 
	  argindex[i] > modify->compute[n]->size_per_molecule_cols)
	error->all(FLERR,
		   "Dump molecule compute vector is accessed out-of-range");

      field2index[i] = add_compute(suffix);
      delete [] suffix;
      
    // fix value = f_ID
    // if no trailing [], then arg is set to 0, else arg is between []

    } else if (strncmp(arg[iarg],"f_",2) == 0) {
      pack_choice[i] = &DumpMolecule::pack_fix;
      vtype[i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all(FLERR,"Invalid attribute in dump molecule command");
	argindex[i] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[i] = 0;

      n = modify->find_fix(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump molecule fix ID");
      if (modify->fix[n]->per_molecule_flag == 0)
	error->all(FLERR,"Dump molecule fix does not compute "
		   "per-molecule info");
      if (argindex[i] == 0 && modify->fix[n]->size_per_molecule_cols > 0)
	error->all(FLERR,"Dump molecule fix does not compute "
		   "per-molecule vector");
      if (argindex[i] > 0 && modify->fix[n]->size_per_molecule_cols == 0)
	error->all(FLERR,"Dump molecule fix does not compute "
		   "per-molecule array");
      if (argindex[i] > 0 && 
	  argindex[i] > modify->fix[n]->size_per_molecule_cols)
	error->all(FLERR,"Dump molecule fix vector is accessed out-of-range");

      field2index[i] = add_fix(suffix);
      delete [] suffix;

    // variable value = v_name

    } else if (strncmp(arg[iarg],"v_",2) == 0) {
      pack_choice[i] = &DumpMolecule::pack_variable;
      vtype[i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      argindex[i] = 0;

      n = input->variable->find(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump molecule variable name");
      if (input->variable->molecule_style(n) == 0)
	error->all(FLERR,"Dump molecule variable is not "
		   "molecule-style variable");

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

int DumpMolecule::add_compute(char *id)
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

int DumpMolecule::add_fix(char *id)
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

int DumpMolecule::add_variable(char *id)
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

/* ---------------------------------------------------------------------- */

int DumpMolecule::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"thresh") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) {
      if (nthresh) {
	memory->destroy(thresh_array);
	memory->destroy(thresh_op);
	memory->destroy(thresh_value);
	thresh_array = NULL;
	thresh_op = NULL;
	thresh_value = NULL;
      }
      nthresh = 0;
      return 2;
    }
    
    if (narg < 4) error->all(FLERR,"Illegal dump_modify command");
    
    // grow threshhold arrays
    
    memory->grow(thresh_array,nthresh+1,"dump:thresh_array");
    memory->grow(thresh_op,(nthresh+1),"dump:thresh_op");
    memory->grow(thresh_value,(nthresh+1),"dump:thresh_value");

    // set attribute type of threshhold
    // customize by adding to if statement
    
    if (strcmp(arg[1],"id") == 0) thresh_array[nthresh] = ID;
    else if (strcmp(arg[1],"type") == 0) thresh_array[nthresh] = TYPE;

    else if (strcmp(arg[1],"x") == 0) thresh_array[nthresh] = X;
    else if (strcmp(arg[1],"y") == 0) thresh_array[nthresh] = Y;
    else if (strcmp(arg[1],"z") == 0) thresh_array[nthresh] = Z;

    else if (strcmp(arg[1],"xs") == 0) thresh_array[nthresh] = XS;
    else if (strcmp(arg[1],"ys") == 0) thresh_array[nthresh] = YS;
    else if (strcmp(arg[1],"zs") == 0) thresh_array[nthresh] = ZS;

    else if (strcmp(arg[1],"vx") == 0) thresh_array[nthresh] = VX;
    else if (strcmp(arg[1],"vy") == 0) thresh_array[nthresh] = VY;
    else if (strcmp(arg[1],"vz") == 0) thresh_array[nthresh] = VZ;

    // compute value = c_ID
    // if no trailing [], then arg is set to 0, else arg is between []
    // must grow field2index and argindex arrays, since access is beyond nfield

    else if (strncmp(arg[1],"c_",2) == 0) {
      thresh_array[nthresh] = COMPUTE;
      memory->grow(field2index,nfield+nthresh+1,"dump:field2index");
      memory->grow(argindex,nfield+nthresh+1,"dump:argindex");
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);
    
      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all(FLERR,"Invalid attribute in dump modify command");
	argindex[nfield+nthresh] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[nfield+nthresh] = 0;
      
      n = modify->find_compute(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump modify compute ID");

      if (modify->compute[n]->per_molecule_flag == 0)
	error->all(FLERR,
		   "Dump modify compute ID does not compute per-molecule info");
      if (argindex[nfield+nthresh] == 0 && 
	  modify->compute[n]->size_per_molecule_cols > 0)
	error->all(FLERR,
		   "Dump modify compute ID does not compute "
		   "per-molecule vector");
      if (argindex[nfield+nthresh] > 0 && 
	  modify->compute[n]->size_per_molecule_cols == 0)
	error->all(FLERR,
		   "Dump modify compute ID does not compute "
		   "per-molecule array");
      if (argindex[nfield+nthresh] > 0 && 
	  argindex[nfield+nthresh] > modify->compute[n]->size_per_molecule_cols)
	error->all(FLERR,"Dump modify compute ID vector is not large enough");

      field2index[nfield+nthresh] = add_compute(suffix);
      delete [] suffix;

    // fix value = f_ID
    // if no trailing [], then arg is set to 0, else arg is between []
    // must grow field2index and argindex arrays, since access is beyond nfield

    } else if (strncmp(arg[1],"f_",2) == 0) {
      thresh_array[nthresh] = FIX;
      memory->grow(field2index,nfield+nthresh+1,"dump:field2index");
      memory->grow(argindex,nfield+nthresh+1,"dump:argindex");
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);
    
      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all(FLERR,"Invalid attribute in dump modify command");
	argindex[nfield+nthresh] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[nfield+nthresh] = 0;
      
      n = modify->find_fix(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump modify fix ID");

      if (modify->fix[n]->per_molecule_flag == 0)
	error->all(FLERR,"Dump modify fix ID does not compute "
		   "per-molecule info");
      if (argindex[nfield+nthresh] == 0 && 
	  modify->fix[n]->size_per_molecule_cols > 0)
	error->all(FLERR,"Dump modify fix ID does not compute "
		   "per-molecule vector");
      if (argindex[nfield+nthresh] > 0 && 
	  modify->fix[n]->size_per_molecule_cols == 0)
	error->all(FLERR,"Dump modify fix ID does not compute "
		   "per-molecule array");
      if (argindex[nfield+nthresh] > 0 && 
	  argindex[nfield+nthresh] > modify->fix[n]->size_per_molecule_cols)
	error->all(FLERR,"Dump modify fix ID vector is not large enough");

      field2index[nfield+nthresh] = add_fix(suffix);
      delete [] suffix;

    // variable value = v_ID
    // must grow field2index and argindex arrays, since access is beyond nfield

    } else if (strncmp(arg[1],"v_",2) == 0) {
      thresh_array[nthresh] = VARIABLE;
      memory->grow(field2index,nfield+nthresh+1,"dump:field2index");
      memory->grow(argindex,nfield+nthresh+1,"dump:argindex");
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);
    
      argindex[nfield+nthresh] = 0;
      
      n = input->variable->find(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump modify variable name");
      if (input->variable->molecule_style(n) == 0)
	error->all(FLERR,"Dump modify variable is not molecule-style variable");

      field2index[nfield+nthresh] = add_variable(suffix);
      delete [] suffix;

    } else error->all(FLERR,"Invalid dump_modify threshhold operator");

    // set operation type of threshhold

    if (strcmp(arg[2],"<") == 0) thresh_op[nthresh] = LT;
    else if (strcmp(arg[2],"<=") == 0) thresh_op[nthresh] = LE;
    else if (strcmp(arg[2],">") == 0) thresh_op[nthresh] = GT;
    else if (strcmp(arg[2],">=") == 0) thresh_op[nthresh] = GE;
    else if (strcmp(arg[2],"==") == 0) thresh_op[nthresh] = EQ;
    else if (strcmp(arg[2],"!=") == 0) thresh_op[nthresh] = NEQ;
    else error->all(FLERR,"Invalid dump_modify threshhold operator");

    // set threshhold value

    thresh_value[nthresh] = atof(arg[3]);

    nthresh++;
    return 4;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf, choose, variable arrays
------------------------------------------------------------------------- */

bigint DumpMolecule::memory_usage()
{
  bigint bytes = Dump::memory_usage();
  bytes += memory->usage(choose,maxlocal);
  bytes += memory->usage(dchoose,maxlocal);
  bytes += memory->usage(clist,maxlocal);
  return bytes;
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
------------------------------------------------------------------------- */

void DumpMolecule::pack_compute(int n)
{
  double *vector = compute[field2index[n]]->vector_molecule;
  double **array = compute[field2index[n]]->array_molecule;
  int index = argindex[n];

  if (index == 0) {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = vector[clist[i]];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = array[clist[i]][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::pack_fix(int n)
{
  double *vector = fix[field2index[n]]->vector_molecule;
  double **array = fix[field2index[n]]->array_molecule;
  int index = argindex[n];

  if (index == 0) {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = vector[clist[i]];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = array[clist[i]][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::pack_variable(int n)
{
  double *vector = vbuf[field2index[n]];

  for (int i = 0; i < nchoose; i++) {
    buf[n] = vector[clist[i]];
    n += size_one;
  }
}

/* ----------------------------------------------------------------------
   one method for every attribute dump molecule can output
   the molecule property is packed into buf starting at n with stride size_one
   customize a new attribute by adding a method
------------------------------------------------------------------------- */

void DumpMolecule::pack_id(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].id;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::pack_type(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].ispecies + 1;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::pack_x(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].x[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::pack_y(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].x[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::pack_z(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].x[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::pack_xs(int n)
{
  Particle::OnePart *particles = particle->particles;

  double boxxlo = domain->boxlo[0];
  double invxprd = 1.0/domain->xprd;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (particles[clist[i]].x[0] - boxxlo) * invxprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::pack_ys(int n)
{
  Particle::OnePart *particles = particle->particles;

  double boxylo = domain->boxlo[1];
  double invyprd = 1.0/domain->yprd;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (particles[clist[i]].x[1] - boxylo) * invyprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::pack_zs(int n)
{
  Particle::OnePart *particles = particle->particles;

  double boxzlo = domain->boxlo[2];
  double invzprd = 1.0/domain->zprd;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (particles[clist[i]].x[2] - boxzlo) * invzprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::pack_vx(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].v[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::pack_vy(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].v[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolecule::pack_vz(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].v[2];
    n += size_one;
  }
}
