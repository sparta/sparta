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
#include "dump_particle.h"
#include "update.h"
#include "domain.h"
#include "particle.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

// customize by adding keyword

enum{ID,TYPE,X,Y,Z,XS,YS,ZS,VX,VY,VZ};
enum{LT,LE,GT,GE,EQ,NEQ};
enum{INT,DOUBLE,STRING};

enum{PERIODIC,OUTFLOW,SPECULAR};            // same as in Domain

/* ---------------------------------------------------------------------- */

DumpParticle::DumpParticle(DSMC *dsmc, int narg, char **arg) :
  Dump(dsmc, narg, arg)
{
  if (narg == 4) error->all(FLERR,"No dump custom arguments specified");

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

  // process attributes
  // ioptional = start of additional optional args
  // only dump image style processes optional args

  ioptional = parse_fields(narg,arg);

  if (ioptional < narg && strcmp(style,"image") != 0)
    error->all(FLERR,"Invalid attribute in dump custom command");
  size_one = nfield = ioptional - 4;

  // atom selection arrays

  maxlocal = 0;
  choose = NULL;
  dchoose = NULL;
  clist = NULL;

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

DumpParticle::~DumpParticle()
{
  delete [] pack_choice;
  delete [] vtype;
  memory->destroy(field2index);
  memory->destroy(argindex);

  memory->destroy(thresh_array);
  memory->destroy(thresh_op);
  memory->destroy(thresh_value);

  memory->destroy(choose);
  memory->destroy(dchoose);
  memory->destroy(clist);

  for (int i = 0; i < size_one; i++) delete [] vformat[i];
  delete [] vformat;

  delete [] columns;
}

/* ---------------------------------------------------------------------- */

void DumpParticle::init_style()
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

  if (binary) header_choice = &DumpParticle::header_binary;
  else header_choice = &DumpParticle::header_item;

  if (binary) write_choice = &DumpParticle::write_binary;
  else write_choice = &DumpParticle::write_text;

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpParticle::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpParticle::header_binary(bigint ndump)
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

void DumpParticle::header_item(bigint ndump)
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

int DumpParticle::count()
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
  }

  // choose all local atoms for output

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
      }

      // unselect atoms that don't meet threshhold criterion

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
  // nchoose = # of selected atoms
  // clist[i] = local index of each selected atom

  nchoose = 0;
  for (i = 0; i < nlocal; i++)
    if (choose[i]) clist[nchoose++] = i;

  return nchoose;
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack()
{
  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);
}

/* ---------------------------------------------------------------------- */

void DumpParticle::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpParticle::write_binary(int n, double *mybuf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(mybuf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpParticle::write_text(int n, double *mybuf)
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

int DumpParticle::parse_fields(int narg, char **arg)
{
  // customize by adding to if statement

  int i;
  for (int iarg = 4; iarg < narg; iarg++) {
    i = iarg-4;

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[i] = &DumpParticle::pack_id;
      vtype[i] = INT;
    } else if (strcmp(arg[iarg],"type") == 0) {
      pack_choice[i] = &DumpParticle::pack_type;
      vtype[i] = INT;

    } else if (strcmp(arg[iarg],"x") == 0) {
      pack_choice[i] = &DumpParticle::pack_x;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"y") == 0) {
      pack_choice[i] = &DumpParticle::pack_y;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"z") == 0) {
      pack_choice[i] = &DumpParticle::pack_z;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xs") == 0) {
      pack_choice[i] = &DumpParticle::pack_xs;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"ys") == 0) {
      pack_choice[i] = &DumpParticle::pack_ys;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"zs") == 0) {
      pack_choice[i] = &DumpParticle::pack_zs;
      vtype[i] = DOUBLE;

    } else if (strcmp(arg[iarg],"vx") == 0) {
      pack_choice[i] = &DumpParticle::pack_vx;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      pack_choice[i] = &DumpParticle::pack_vy;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      pack_choice[i] = &DumpParticle::pack_vz;
      vtype[i] = DOUBLE;

    } else return iarg;
  }

  return narg;
}

/* ---------------------------------------------------------------------- */

int DumpParticle::modify_param(int narg, char **arg)
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

    else error->all(FLERR,"Invalid dump_modify threshhold operator");

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

bigint DumpParticle::memory_usage()
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

/* ----------------------------------------------------------------------
   one method for every attribute dump custom can output
   the atom property is packed into buf starting at n with stride size_one
   customize a new attribute by adding a method
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_id(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].id;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_type(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].ispecies + 1;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_x(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].x[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_y(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].x[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_z(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].x[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_xs(int n)
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

void DumpParticle::pack_ys(int n)
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

void DumpParticle::pack_zs(int n)
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

void DumpParticle::pack_vx(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].v[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_vy(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].v[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_vz(int n)
{
  Particle::OnePart *particles = particle->particles;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = particles[clist[i]].v[2];
    n += size_one;
  }
}
