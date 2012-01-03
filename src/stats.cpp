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

#include "dsmctype.h"
#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "stats.h"
#include "update.h"
#include "particle.h"
#include "domain.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

// customize a new keyword by adding to this list:

// step, elapsed, dt, cpu, tpcpu, spcpu
// npart, temp, vol, lx, ly, lz, xlo, xhi, ylo, yhi, zlo, zhi

enum{INT,FLOAT,BIGINT};

#define MAXLINE 8192               // make this 4x longer than Input::MAXLINE

/* ---------------------------------------------------------------------- */

Stats::Stats(DSMC *dsmc) : Pointers(dsmc)
{
  MPI_Comm_rank(world,&me);

  line = new char[MAXLINE];

  keyword = NULL;
  vfunc = NULL;
  vtype = NULL;

  format = NULL;
  format_user = NULL;

  field2index = NULL;
  argindex1 = NULL;
  argindex2 = NULL;

  // default args

  char **arg = new char*[4];
  arg[0] = (char *) "step";
  arg[1] = (char *) "cpu";
  arg[2] = (char *) "npart";
  arg[3] = (char *) "temp";

  nfield = 0;
  set_fields(4,arg);

  delete [] arg;

  // stats_modify defaults

  flushflag = 0;

  // format strings

  char *bigint_format = (char *) BIGINT_FORMAT;

  format_float_one_def = (char *) "%12.8g";
  format_int_one_def = (char *) "%8d";
  sprintf(format_bigint_one_def,"%%8%s",&bigint_format[1]);

  format_float_user = NULL;
  format_int_user = NULL;
  format_bigint_user = NULL;
}

/* ---------------------------------------------------------------------- */

Stats::~Stats()
{
  delete [] line;
  deallocate();

  // format strings

  delete [] format_float_user;
  delete [] format_int_user;
  delete [] format_bigint_user;
}

/* ---------------------------------------------------------------------- */

void Stats::init()
{
  // set format string for each field
  // add trailing '/n' to last value

  char *ptr;
  for (int i = 0; i < nfield; i++) {
    format[i][0] = '\0';

    if (format_user[i]) ptr = format_user[i];
    else if (vtype[i] == FLOAT) {
      if (format_float_user) ptr = format_float_user;
      else ptr = format_float_one_def;
    } else if (vtype[i] == INT) {
      if (format_int_user) ptr = format_int_user;
      else ptr = format_int_one_def;
    } else if (vtype[i] == BIGINT) {
      if (format_bigint_user) ptr = format_bigint_user;
      else ptr = format_bigint_one_def;
    }

    int n = strlen(format[i]);
    sprintf(&format[i][n],"%s ",ptr);

    if (i == nfield-1) strcat(format[i],"\n");
  }
}

/* ---------------------------------------------------------------------- */

void Stats::header()
{
  int loc = 0;
  for (int i = 0; i < nfield; i++)
    loc += sprintf(&line[loc],"%s ",keyword[i]);
  sprintf(&line[loc],"\n");
  
  if (me == 0) {
    if (screen) fprintf(screen,line);
    if (logfile) fprintf(logfile,line);
  }
}

/* ---------------------------------------------------------------------- */

void Stats::compute(int flag)
{
  int i;

  firststep = flag;
  bigint ntimestep = update->ntimestep;

  // add each stat value to line with its specific format

  int loc = 0;
  for (ifield = 0; ifield < nfield; ifield++) {
    (this->*vfunc[ifield])();
    if (vtype[ifield] == FLOAT)
      loc += sprintf(&line[loc],format[ifield],dvalue);
    else if (vtype[ifield] == INT) 
      loc += sprintf(&line[loc],format[ifield],ivalue);
    else if (vtype[ifield] == BIGINT) {
      loc += sprintf(&line[loc],format[ifield],bivalue);
    }
  }

  // print line to screen and logfile

  if (me == 0) {
    if (screen) fprintf(screen,line);
    if (logfile) {
      fprintf(logfile,line);
      if (flushflag) fflush(logfile);
    }
  }
}

/* ----------------------------------------------------------------------
   modify stats parameters
------------------------------------------------------------------------- */

void Stats::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal stats_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal stats_modify command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
	delete [] output->var_stats;
	int n = strlen(&arg[iarg+1][2]) + 1;
	output->var_stats = new char[n];
	strcpy(output->var_stats,&arg[iarg+1][2]);
      } else error->all(FLERR,"Illegal stats_modify command");
      output->stats_every = 0;
      iarg += 2;

    } else if (strcmp(arg[iarg],"flush") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal stats_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) flushflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) flushflag = 1;
      else error->all(FLERR,"Illegal stats_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal stats_modify command");
      if (strcmp(arg[iarg+1],"int") == 0) {
	if (format_int_user) delete [] format_int_user;
	int n = strlen(arg[iarg+2]) + 1;
	format_int_user = new char[n];
	strcpy(format_int_user,arg[iarg+2]);
	if (format_bigint_user) delete [] format_bigint_user;
	n = strlen(format_int_user) + 3;
	format_bigint_user = new char[n];
	char *ptr = strchr(format_int_user,'d');
	if (ptr == NULL) 
	  error->all(FLERR,
		     "Stats_modify int format does not contain d character");
	*ptr = '\0';
	sprintf(format_bigint_user,"%s%s%s",format_int_user,
		BIGINT_FORMAT,ptr+1);
	*ptr = 'd';
      } else if (strcmp(arg[iarg+1],"float") == 0) {
	if (format_float_user) delete [] format_float_user;
	int n = strlen(arg[iarg+2]) + 1;
	format_float_user = new char[n];
	strcpy(format_float_user,arg[iarg+2]);
      } else {
	int i = atoi(arg[iarg+1]) - 1;
	if (i < 0 || i >= nfield)
	  error->all(FLERR,"Illegal stats_modify command");
	if (format_user[i]) delete [] format_user[i];
	int n = strlen(arg[iarg+2]) + 1;
	format_user[i] = new char[n];
	strcpy(format_user[i],arg[iarg+2]);
      }
      iarg += 3;

    } else error->all(FLERR,"Illegal stats_modify command");
  }
}

/* ----------------------------------------------------------------------
   allocate all per-field memory
------------------------------------------------------------------------- */

void Stats::allocate()
{
  int n = nfield;

  keyword = new char*[n];
  for (int i = 0; i < n; i++) keyword[i] = new char[32];
  vfunc = new FnPtr[n];
  vtype = new int[n];

  format = new char*[n];
  for (int i = 0; i < n; i++) format[i] = new char[32];
  format_user = new char*[n];
  for (int i = 0; i < n; i++) format_user[i] = NULL;

  field2index = new int[n];
  argindex1 = new int[n];
  argindex2 = new int[n];
}

/* ----------------------------------------------------------------------
   deallocate all per-field memory
------------------------------------------------------------------------- */

void Stats::deallocate()
{
  int n = nfield;

  for (int i = 0; i < n; i++) delete [] keyword[i];
  delete [] keyword;
  delete [] vfunc;
  delete [] vtype;

  for (int i = 0; i < n; i++) delete [] format[i];
  delete [] format;
  for (int i = 0; i < n; i++) delete [] format_user[i];
  delete [] format_user;

  delete [] field2index;
  delete [] argindex1;
  delete [] argindex2;
}

/* ----------------------------------------------------------------------
   set fields of stats output from args
------------------------------------------------------------------------- */

void Stats::set_fields(int narg, char **arg)
{
  deallocate();
  nfield = narg;
  allocate();

  nfield = 0;

  // customize a new keyword by adding to if statement

  for (int i = 0; i < narg; i++) {
    if (strcmp(arg[i],"step") == 0) {
      addfield("Step",&Stats::compute_step,BIGINT);
    } else if (strcmp(arg[i],"elapsed") == 0) {
      addfield("Elapsed",&Stats::compute_elapsed,BIGINT);
    } else if (strcmp(arg[i],"dt") == 0) {
      addfield("Dt",&Stats::compute_dt,FLOAT);
    } else if (strcmp(arg[i],"cpu") == 0) {
      addfield("CPU",&Stats::compute_cpu,FLOAT);
    } else if (strcmp(arg[i],"tpcpu") == 0) {
      addfield("T/CPU",&Stats::compute_tpcpu,FLOAT);
    } else if (strcmp(arg[i],"spcpu") == 0) {
      addfield("S/CPU",&Stats::compute_spcpu,FLOAT);

    } else if (strcmp(arg[i],"npart") == 0) {
      addfield("Npart",&Stats::compute_npart,BIGINT);
    } else if (strcmp(arg[i],"temp") == 0) {
      addfield("Temp",&Stats::compute_temp,FLOAT);

    } else if (strcmp(arg[i],"vol") == 0) {
      addfield("Volume",&Stats::compute_vol,FLOAT);
    } else if (strcmp(arg[i],"lx") == 0) {
      addfield("Lx",&Stats::compute_lx,FLOAT);
    } else if (strcmp(arg[i],"ly") == 0) {
      addfield("Ly",&Stats::compute_ly,FLOAT);
    } else if (strcmp(arg[i],"lz") == 0) {
      addfield("Lz",&Stats::compute_lz,FLOAT);

    } else if (strcmp(arg[i],"xlo") == 0) {
      addfield("Xlo",&Stats::compute_xlo,FLOAT);
    } else if (strcmp(arg[i],"xhi") == 0) {
      addfield("Xhi",&Stats::compute_xhi,FLOAT);
    } else if (strcmp(arg[i],"ylo") == 0) {
      addfield("Ylo",&Stats::compute_ylo,FLOAT);
    } else if (strcmp(arg[i],"yhi") == 0) {
      addfield("Yhi",&Stats::compute_yhi,FLOAT);
    } else if (strcmp(arg[i],"zlo") == 0) {
      addfield("Zlo",&Stats::compute_zlo,FLOAT);
    } else if (strcmp(arg[i],"zhi") == 0) {
      addfield("Zhi",&Stats::compute_zhi,FLOAT);

    } else error->all(FLERR,"Invalid keyword in stats_style command");
  }
}

/* ----------------------------------------------------------------------
   add field to list of quantities to print
------------------------------------------------------------------------- */

void Stats::addfield(const char *key, FnPtr func, int typeflag)
{
  strcpy(keyword[nfield],key);
  vfunc[nfield] = func;
  vtype[nfield] = typeflag;
  nfield++;
}

/* ----------------------------------------------------------------------
   compute a single thermodyanmic value, word is any keyword in custom list
   called when a variable is evaluated by Variable class
   return value as double in answer
   return 0 if str is recoginzed keyword, 1 if unrecognized
   customize a new keyword by adding to if statement
------------------------------------------------------------------------- */

int Stats::evaluate_keyword(char *word, double *answer)
{
  return 0;
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
   compute/fix are normalized by atoms if returning extensive value
   variable value is not normalized (formula should normalize if desired)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   one method for every keyword thermo can output
   called by compute() or evaluate_keyword()
   compute will have already been called
   set ivalue/dvalue/bivalue if value is int/double/bigint
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

void Stats::compute_step()
{
  bivalue = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_elapsed()
{
  bivalue = update->ntimestep - update->firststep;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_dt()
{
  dvalue = update->dt;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_cpu()
{
  if (firststep == 0) dvalue = 0.0;
  else dvalue = timer->elapsed(TIME_LOOP);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_tpcpu()
{
  double new_cpu;
  double new_time = update->ntimestep * update->dt;

  if (firststep == 0) {
    new_cpu = 0.0;
    dvalue = 0.0;
  } else {
    new_cpu = timer->elapsed(TIME_LOOP);
    double cpu_diff = new_cpu - last_tpcpu;
    double time_diff = new_time - last_time;
    if (time_diff > 0.0 && cpu_diff > 0.0) dvalue = time_diff/cpu_diff;
    else dvalue = 0.0;
  }

  last_time = new_time;
  last_tpcpu = new_cpu;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_spcpu()
{
  double new_cpu;
  int new_step = update->ntimestep;

  if (firststep == 0) {
    new_cpu = 0.0;
    dvalue = 0.0;
  } else {
    new_cpu = timer->elapsed(TIME_LOOP);
    double cpu_diff = new_cpu - last_spcpu;
    int step_diff = new_step - last_step;
    if (cpu_diff > 0.0) dvalue = step_diff/cpu_diff;
    else dvalue = 0.0;
  }

  last_step = new_step;
  last_spcpu = new_cpu;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_npart()
{
  bigint n = particle->nlocal;
  MPI_Allreduce(&n,&particle->nglobal,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  bivalue = particle->nglobal;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_temp()
{
  dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_vol()
{
  if (domain->dimension == 3)
    dvalue = domain->xprd * domain->yprd * domain->zprd;
  else
    dvalue = domain->xprd * domain->yprd;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_lx()
{
  dvalue = domain->xprd;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_ly()
{
  dvalue = domain->yprd;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_lz()
{
  dvalue = domain->zprd;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_xlo()
{
  dvalue = domain->boxlo[0];
}

/* ---------------------------------------------------------------------- */

void Stats::compute_xhi()
{
  dvalue = domain->boxhi[0];
}

/* ---------------------------------------------------------------------- */

void Stats::compute_ylo()
{
  dvalue = domain->boxlo[1];
}

/* ---------------------------------------------------------------------- */

void Stats::compute_yhi()
{
  dvalue = domain->boxhi[1];
}

/* ---------------------------------------------------------------------- */

void Stats::compute_zlo()
{
  dvalue = domain->boxlo[2];
}

/* ---------------------------------------------------------------------- */

void Stats::compute_zhi()
{
  dvalue = domain->boxhi[2];
}

