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

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "output.h"
#include "style_dump.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "particle.h"
#include "grid.h"
#include "surf.h"
#include "domain.h"
#include "modify.h"
#include "stats.h"
#include "dump.h"
#include "write_restart.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

#define DELTA 1

/* ----------------------------------------------------------------------
   initialize all output 
------------------------------------------------------------------------- */

Output::Output(DSMC *dsmc) : Pointers(dsmc)
{
  // create default Stats class

  stats = new Stats(dsmc);
    
  stats_every = 0;
  var_stats = NULL;
  
  ndump = 0;
  max_dump = 0;
  every_dump = NULL;
  next_dump = NULL;
  last_dump = NULL;
  var_dump = NULL;
  ivar_dump = NULL;
  dump = NULL;

  restart = NULL;
  restart1 = restart2 = NULL;
  restart_every = 0;
  last_restart = -1;
}

/* ----------------------------------------------------------------------
   free all memory 
------------------------------------------------------------------------- */

Output::~Output()
{
  if (stats) delete stats;
  delete [] var_stats;

  memory->destroy(every_dump);
  memory->destroy(next_dump);
  memory->destroy(last_dump);
  for (int i = 0; i < ndump; i++) delete [] var_dump[i];
  memory->sfree(var_dump);
  memory->destroy(ivar_dump);
  for (int i = 0; i < ndump; i++) delete dump[i];
  memory->sfree(dump);

  delete restart;
  delete [] restart1;
  delete [] restart2;
}

/* ---------------------------------------------------------------------- */

void Output::init()
{
  stats->init();
  if (stats_every) delete [] var_stats;
  else if (var_stats) {
    ivar_stats = input->variable->find(var_stats);
    if (ivar_stats < 0)
      error->all(FLERR,"Variable name for stats every does not exist");
    if (!input->variable->equal_style(ivar_stats))
      error->all(FLERR,"Variable for stats every is invalid style");
  }

  for (int i = 0; i < ndump; i++) dump[i]->init();
  for (int i = 0; i < ndump; i++)
    if (every_dump[i] == 0) {
      ivar_dump[i] = input->variable->find(var_dump[i]);
      if (ivar_dump[i] < 0)
	error->all(FLERR,"Variable name for dump every does not exist");
      if (!input->variable->equal_style(ivar_dump[i]))
	error->all(FLERR,"Variable for dump every is invalid style");
    }
}

/* ----------------------------------------------------------------------
   perform output for setup of run/min
   do dump first, so memory_usage will include dump allocation
   do stats last, so will print after memory_usage
------------------------------------------------------------------------- */

void Output::setup(int flag)
{
  bigint ntimestep = update->ntimestep;

  // perform dump at start of run if current timestep is multiple of every
  //   and last dump was not on this timestep
  // set next_dump to multiple of every
  // will not write on last step of run unless multiple of every
  // set next_dump_any to smallest next_dump
  // if no dumps, set next_dump_any to last+1 so will not influence next
  // wrap dumps that invoke computes with clear/add
  // if dump not written now, add_all on future step since clear/add is noop

  int writeflag;

  if (ndump) {
    for (int idump = 0; idump < ndump; idump++) {
      if (dump[idump]->clearstep) modify->clearstep_compute();
      writeflag = 0;
      if (every_dump[idump] && ntimestep % every_dump[idump] == 0 && 
	  last_dump[idump] != ntimestep) writeflag = 1;
      if (last_dump[idump] < 0 && dump[idump]->first_flag == 1) writeflag = 1;
      if (writeflag) {
	dump[idump]->write();
	last_dump[idump] = ntimestep;
      }
      if (every_dump[idump])
	next_dump[idump] = 
	  (ntimestep/every_dump[idump])*every_dump[idump] + every_dump[idump];
      else {
	int nextdump = static_cast<int> 
	  (input->variable->compute_equal(ivar_dump[idump]));
	if (nextdump <= ntimestep)
	  error->all(FLERR,"Dump every variable returned a bad timestep");
	next_dump[idump] = nextdump;
      }
      if (dump[idump]->clearstep) {
	if (writeflag) modify->addstep_compute(next_dump[idump]);
	else modify->addstep_compute_all(next_dump[idump]);
      }
      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  } else next_dump_any = update->laststep + 1;

  // do not write a restart file at start of run
  // set next_restart to multiple of every
  // will not write on last step of run unless multiple of every
  // if every = 0, set next_restart to last+1 so will not influence next

  if (restart_every)
    next_restart = (ntimestep/restart_every)*restart_every + restart_every;
  else next_restart = update->laststep + 1;

  // print memory usage unless being called between multiple runs

  if (flag) memory_usage();

  // always do stats with header at start of run
  // set next_stats to multiple of every or last step of run (if smaller)
  // if every = 0, set next_stats to last step of run
  // stats may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  stats->header();
  stats->compute(0);
  last_stats = ntimestep;

  if (stats_every) {
    next_stats = (ntimestep/stats_every)*stats_every + stats_every;
    next_stats = MIN(next_stats,update->laststep);
  } else if (var_stats) {
    next_stats = static_cast<int> 
      (input->variable->compute_equal(ivar_stats));
    if (next_stats <= ntimestep)
      error->all(FLERR,"stats every variable returned a bad timestep");
  } else next_stats = update->laststep;

  //modify->addstep_compute(next_stats);

  // next = next timestep any output will be done

  next = MIN(next_dump_any,next_restart);
  next = MIN(next,next_stats);
}

/* ----------------------------------------------------------------------
   perform all output for this timestep
   only perform output if next matches current step and last doesn't
   do dump/restart before stats so stats CPU time will include them
------------------------------------------------------------------------- */

void Output::write(bigint ntimestep)
{
  // next_dump does not force output on last step of run
  // wrap dumps that invoke computes with clear/add
  // download data from GPU if necessary

  if (next_dump_any == ntimestep) {

    for (int idump = 0; idump < ndump; idump++) {
      if (next_dump[idump] == ntimestep && last_dump[idump] != ntimestep) {
        if (dump[idump]->clearstep) modify->clearstep_compute();
	dump[idump]->write();
	last_dump[idump] = ntimestep;
	if (every_dump[idump]) next_dump[idump] += every_dump[idump];
	else {
	  int nextdump = static_cast<int> 
	    (input->variable->compute_equal(ivar_dump[idump]));
	  if (nextdump <= ntimestep)
	    error->all(FLERR,"Dump every variable returned a bad timestep");
	  next_dump[idump] = nextdump;
	}
        if (dump[idump]->clearstep) modify->addstep_compute(next_dump[idump]);
      }
      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  }

  // next_restart does not force output on last step of run
  // for toggle = 0, replace "*" with current timestep in restart filename
  // download data from GPU if necessary

  if (next_restart == ntimestep && last_restart != ntimestep) {

    if (restart_toggle == 0) {
      char *file = new char[strlen(restart1) + 16];
      char *ptr = strchr(restart1,'*');
      *ptr = '\0';
      sprintf(file,"%s" BIGINT_FORMAT "%s",restart1,ntimestep,ptr+1);
      *ptr = '*';
      restart->write(file);
      delete [] file;
    } else if (restart_toggle == 1) {
      restart->write(restart1);
      restart_toggle = 2;
    } else if (restart_toggle == 2) {
      restart->write(restart2);
      restart_toggle = 1;
    }
    last_restart = ntimestep;
    next_restart += restart_every;
  }

  // insure next_stats forces output on last step of run
  // stats may invoke computes so wrap with clear/add

  if (next_stats == ntimestep && last_stats != ntimestep) {
    modify->clearstep_compute();
    stats->compute(1);
    last_stats = ntimestep;
    if (stats_every) next_stats += stats_every;
    else if (var_stats) {
      next_stats = static_cast<int> 
	(input->variable->compute_equal(ivar_stats));
      if (next_stats <= ntimestep)
	error->all(FLERR,"stats every variable returned a bad timestep");
    } else next_stats = update->laststep;
    next_stats = MIN(next_stats,update->laststep);
    modify->addstep_compute(next_stats);
  }

  // next = next timestep any output will be done

  next = MIN(next_dump_any,next_restart);
  next = MIN(next,next_stats);
}

/* ----------------------------------------------------------------------
   force a snapshot to be written for all dumps
------------------------------------------------------------------------- */

void Output::write_dump(bigint ntimestep)
{
  for (int idump = 0; idump < ndump; idump++) {
    dump[idump]->write();
    last_dump[idump] = ntimestep;
  }
}

/* ----------------------------------------------------------------------
   force a restart file to be written
------------------------------------------------------------------------- */

void Output::write_restart(bigint ntimestep)
{
  if (restart_toggle == 0) {
    char *file = new char[strlen(restart1) + 16];
    char *ptr = strchr(restart1,'*');
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s",restart1,ntimestep,ptr+1);
    *ptr = '*';
    restart->write(file);
    delete [] file;
  } else if (restart_toggle == 1) {
    restart->write(restart1);
    restart_toggle = 2;
  } else if (restart_toggle == 2) {
    restart->write(restart2);
    restart_toggle = 1;
  }

  last_restart = ntimestep;
}

/* ----------------------------------------------------------------------
   add a Dump to list of Dumps 
------------------------------------------------------------------------- */

void Output::add_dump(int narg, char **arg)
{
  if (narg < 4) error->all(FLERR,"Illegal dump command");

  // error checks

  for (int idump = 0; idump < ndump; idump++)
    if (strcmp(arg[0],dump[idump]->id) == 0) 
      error->all(FLERR,"Reuse of dump ID");
  if (atoi(arg[2]) <= 0) error->all(FLERR,"Invalid dump frequency");

  // extend Dump list if necessary

  if (ndump == max_dump) {
    max_dump += DELTA;
    dump = (Dump **)
      memory->srealloc(dump,max_dump*sizeof(Dump *),"output:dump");
    memory->grow(every_dump,max_dump,"output:every_dump");
    memory->grow(next_dump,max_dump,"output:next_dump");
    memory->grow(last_dump,max_dump,"output:last_dump");
    var_dump = (char **)
      memory->srealloc(var_dump,max_dump*sizeof(char *),"output:var_dump");
    memory->grow(ivar_dump,max_dump,"output:ivar_dump");
  }

  // create the Dump

  if (0) return;         // dummy line to enable else-if macro expansion

#define DUMP_CLASS
#define DumpStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) dump[ndump] = new Class(dsmc,narg,arg);
#include "style_dump.h"
#undef DUMP_CLASS

  else error->all(FLERR,"Invalid dump style");

  every_dump[ndump] = atoi(arg[2]);
  if (every_dump[ndump] <= 0) error->all(FLERR,"Illegal dump command");
  last_dump[ndump] = -1;
  var_dump[ndump] = NULL;
  ndump++;
}

/* ----------------------------------------------------------------------
   modify parameters of a Dump 
------------------------------------------------------------------------- */

void Output::modify_dump(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal dump_modify command");

  // find which dump it is

  int idump;
  for (idump = 0; idump < ndump; idump++)
    if (strcmp(arg[0],dump[idump]->id) == 0) break;
  if (idump == ndump) error->all(FLERR,"Cound not find dump_modify ID");

  dump[idump]->modify_params(narg-1,&arg[1]);
}

/* ----------------------------------------------------------------------
   delete a Dump from list of Dumps 
------------------------------------------------------------------------- */

void Output::delete_dump(char *id)
{
  // find which dump it is and delete it

  int idump;
  for (idump = 0; idump < ndump; idump++)
    if (strcmp(id,dump[idump]->id) == 0) break;
  if (idump == ndump) error->all(FLERR,"Could not find undump ID");

  delete dump[idump];
  delete [] var_dump[idump];

  // move other dumps down in list one slot

  for (int i = idump+1; i < ndump; i++) {
    dump[i-1] = dump[i];
    every_dump[i-1] = every_dump[i];
    next_dump[i-1] = next_dump[i];
    last_dump[i-1] = last_dump[i];
    var_dump[i-1] = var_dump[i];
    ivar_dump[i-1] = ivar_dump[i];
  }
  ndump--;
}

/* ----------------------------------------------------------------------
   new stats style 
------------------------------------------------------------------------- */

void Output::create_stats(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal stats_style command");
  stats->set_fields(narg,arg);
}

/* ----------------------------------------------------------------------
   setup restart capability
   if only one filename and it contains no "*", then append ".*"
------------------------------------------------------------------------- */

void Output::create_restart(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal restart command");

  if (restart) delete restart;
  delete [] restart1;
  delete [] restart2;
  restart = NULL;
  restart1 = restart2 = NULL;
  last_restart = -1;

  restart_every = atoi(arg[0]);
  if (restart_every == 0) {
    if (narg != 1) error->all(FLERR,"Illegal restart command");
    return;
  }

  restart = new WriteRestart(dsmc);

  int n = strlen(arg[1]) + 3;
  restart1 = new char[n];
  strcpy(restart1,arg[1]);

  if (narg == 2) {
    restart_toggle = 0;
    restart2 = NULL;
    if (strchr(restart1,'*') == NULL) strcat(restart1,".*");
  } else if (narg == 3) {
    restart_toggle = 1;
    n = strlen(arg[2]) + 1;
    restart2 = new char[n];
    strcpy(restart2,arg[2]);
  } else error->all(FLERR,"Illegal restart command");
}

/* ----------------------------------------------------------------------
   sum and print memory usage
------------------------------------------------------------------------- */

void Output::memory_usage()
{
  bigint pbytes,gbytes,sbytes,bytes;
  pbytes = particle->memory_usage();
  gbytes = grid->memory_usage();
  sbytes = surf->memory_usage();
  bytes = pbytes + gbytes + sbytes;
  bytes += modify->memory_usage();

  double scale = 1.0/1024.0/1024.0;

  bigint ave,min,max;

  MPI_Allreduce(&pbytes,&ave,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  double pave = scale * ave/comm->nprocs;
  MPI_Allreduce(&pbytes,&min,1,MPI_DSMC_BIGINT,MPI_MIN,world);
  double pmin = scale * min;
  MPI_Allreduce(&pbytes,&max,1,MPI_DSMC_BIGINT,MPI_MAX,world);
  double pmax = scale * max;

  MPI_Allreduce(&gbytes,&ave,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  double gave = scale * ave/comm->nprocs;
  MPI_Allreduce(&gbytes,&min,1,MPI_DSMC_BIGINT,MPI_MIN,world);
  double gmin = scale * min;
  MPI_Allreduce(&gbytes,&max,1,MPI_DSMC_BIGINT,MPI_MAX,world);
  double gmax = scale * max;

  MPI_Allreduce(&sbytes,&ave,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  double save = scale * ave/comm->nprocs;
  MPI_Allreduce(&sbytes,&min,1,MPI_DSMC_BIGINT,MPI_MIN,world);
  double smin = scale * min;
  MPI_Allreduce(&sbytes,&max,1,MPI_DSMC_BIGINT,MPI_MAX,world);
  double smax = scale * max;

  MPI_Allreduce(&bytes,&ave,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  double tave = scale * ave/comm->nprocs;
  MPI_Allreduce(&bytes,&min,1,MPI_DSMC_BIGINT,MPI_MIN,world);
  double tmin = scale * min;
  MPI_Allreduce(&bytes,&max,1,MPI_DSMC_BIGINT,MPI_MAX,world);
  double tmax = scale * max;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Memory usage per proc in Mbytes:\n");
      fprintf(screen,"  particles (ave,min,max) = %g %g %g\n",
	      pave,pmin,pmax);
      fprintf(screen,"  grid      (ave,min,max) = %g %g %g\n",
	      gave,gmin,gmax);
      fprintf(screen,"  surf      (ave,min,max) = %g %g %g\n",
	      save,smin,smax);
      fprintf(screen,"  total     (ave,min,max) = %g %g %g\n",
	      tave,tmin,tmax);
    }
    if (logfile) {
      fprintf(logfile,"Memory usage per proc in Mbytes:\n");
      fprintf(logfile,"  particles (ave,min,max) = %g %g %g\n",
	      pave,pmin,pmax);
      fprintf(logfile,"  grid      (ave,min,max) = %g %g %g\n",
	      gave,gmin,gmax);
      fprintf(logfile,"  surf      (ave,min,max) = %g %g %g\n",
	      save,smin,smax);
      fprintf(logfile,"  total     (ave,min,max) = %g %g %g\n",
	      tave,tmin,tmax);
    }
  }
}
