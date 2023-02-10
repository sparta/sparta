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

using namespace SPARTA_NS;

#define DELTA 1

/* ----------------------------------------------------------------------
   initialize all output
------------------------------------------------------------------------- */

Output::Output(SPARTA *sparta) : Pointers(sparta)
{
  // create default Stats class

  stats = new Stats(sparta);

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

  restart_flag = restart_flag_single = restart_flag_double = 0;
  restart_every_single = restart_every_double = 0;
  last_restart = -1;
  restart1 = restart2a = restart2b = NULL;
  var_restart_single = var_restart_double = NULL;
  restart = NULL;
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

  delete [] restart1;
  delete [] restart2a;
  delete [] restart2b;
  delete [] var_restart_single;
  delete [] var_restart_double;
  delete restart;
}

/* ---------------------------------------------------------------------- */

void Output::init()
{
  stats->init();
  if (var_stats) {
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

  if (restart_flag_single && restart_every_single == 0) {
    ivar_restart_single = input->variable->find(var_restart_single);
    if (ivar_restart_single < 0)
      error->all(FLERR,"Variable name for restart does not exist");
    if (!input->variable->equal_style(ivar_restart_single))
      error->all(FLERR,"Variable for restart is invalid style");
  }
  if (restart_flag_double && restart_every_double == 0) {
    ivar_restart_double = input->variable->find(var_restart_double);
    if (ivar_restart_double < 0)
      error->all(FLERR,"Variable name for restart does not exist");
    if (!input->variable->equal_style(ivar_restart_double))
      error->all(FLERR,"Variable for restart is invalid style");
  }
}

/* ----------------------------------------------------------------------
   perform output for setup of run/min
   do dump first, so memory_usage will include dump allocation
   do stats last, so will print after memory_usage
   memflag = 0/1 for printing out memory usage
------------------------------------------------------------------------- */

void Output::setup(int memflag)
{
  bigint ntimestep = update->ntimestep;

  // perform dump at start of run only if:
  //   current timestep is multiple of every and last dump not >= this step
  //   this is first run after dump created and firstflag is set
  //   note that variable freq will not write unless triggered by firstflag
  // set next_dump to multiple of every or variable value
  // set next_dump_any to smallest next_dump
  // wrap dumps that invoke computes and variable eval with clear/add
  // if dump not written now, use addstep_compute_all() since don't know
  //   what computes the dump write would invoke
  // if no dumps, set next_dump_any to last+1 so will not influence next

  int writeflag;

  if (ndump) {
    for (int idump = 0; idump < ndump; idump++) {
      if (dump[idump]->clearstep || every_dump[idump] == 0)
        modify->clearstep_compute();
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
        bigint nextdump = static_cast<bigint>
          (input->variable->compute_equal(ivar_dump[idump]));
        if (nextdump <= ntimestep)
          error->all(FLERR,"Dump every variable returned a bad timestep");
        next_dump[idump] = nextdump;
      }
      if (dump[idump]->clearstep || every_dump[idump] == 0) {
        if (writeflag) modify->addstep_compute(next_dump[idump]);
        else modify->addstep_compute_all(next_dump[idump]);
      }
      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  } else next_dump_any = update->laststep + 1;

  // do not write restart files at start of run
  // set next_restart values to multiple of every or variable value
  // wrap variable eval with clear/add
  // if no restarts, set next_restart to last+1 so will not influence next

  if (restart_flag) {
    if (restart_flag_single) {
      if (restart_every_single)
        next_restart_single =
          (ntimestep/restart_every_single)*restart_every_single +
          restart_every_single;
      else {
        bigint nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_single));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad timestep");
        next_restart_single = nextrestart;
      }
    } else next_restart_single = update->laststep + 1;
    if (restart_flag_double) {
      if (restart_every_double)
        next_restart_double =
          (ntimestep/restart_every_double)*restart_every_double +
          restart_every_double;
      else {
        bigint nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_double));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad timestep");
        next_restart_double = nextrestart;
      }
    } else next_restart_double = update->laststep + 1;
    next_restart = MIN(next_restart_single,next_restart_double);
  } else next_restart = update->laststep + 1;

  // print memory usage unless being called between multiple runs

  if (memflag) memory_usage();

  // set next_stats to multiple of every or variable eval if var defined
  // insure stats output on last step of run
  // stats may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  stats->header();
  stats->compute(0);
  last_stats = ntimestep;

  if (var_stats) {
    next_stats = static_cast<bigint>
      (input->variable->compute_equal(ivar_stats));
    if (next_stats <= ntimestep)
      error->all(FLERR,"Stats every variable returned a bad timestep");
  } else if (stats_every) {
    next_stats = (ntimestep/stats_every)*stats_every + stats_every;
    next_stats = MIN(next_stats,update->laststep);
  } else next_stats = update->laststep;

  modify->addstep_compute(next_stats);

  // next = next timestep any output will be done

  next = MIN(next_dump_any,next_restart);
  next = MIN(next,next_stats);
}

/* ----------------------------------------------------------------------
   perform all output for this timestep
   only perform output if next matches current step and last output doesn't
   do dump/restart before stats so stats CPU time will include them
------------------------------------------------------------------------- */

void Output::write(bigint ntimestep)
{
  // next_dump does not force output on last step of run
  // wrap dumps that invoke computes or eval of variable with clear/add
  // download data from GPU if necessary

  if (next_dump_any == ntimestep) {
    for (int idump = 0; idump < ndump; idump++) {
      if (next_dump[idump] == ntimestep) {
        if (dump[idump]->clearstep || every_dump[idump] == 0)
          modify->clearstep_compute();
        if (last_dump[idump] != ntimestep) {
          dump[idump]->write();
          last_dump[idump] = ntimestep;
        }
        if (every_dump[idump]) next_dump[idump] += every_dump[idump];
        else {
          bigint nextdump = static_cast<bigint>
            (input->variable->compute_equal(ivar_dump[idump]));
          if (nextdump <= ntimestep)
            error->all(FLERR,"Dump every variable returned a bad timestep");
          next_dump[idump] = nextdump;
        }
        if (dump[idump]->clearstep || every_dump[idump] == 0)
          modify->addstep_compute(next_dump[idump]);
      }
      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  }

  // next_restart does not force output on last step of run
  // for toggle = 0, replace "*" with current timestep in restart filename
  // eval of variable may invoke computes so wrap with clear/add

  if (next_restart == ntimestep) {
    if (next_restart_single == ntimestep) {
      char *file = new char[strlen(restart1) + 16];
      char *ptr = strchr(restart1,'*');
      *ptr = '\0';
      sprintf(file,"%s" BIGINT_FORMAT "%s",restart1,ntimestep,ptr+1);
      *ptr = '*';
      if (last_restart != ntimestep) restart->write(file);
      delete [] file;
      if (restart_every_single) next_restart_single += restart_every_single;
      else {
        modify->clearstep_compute();
        bigint nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_single));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad timestep");
        next_restart_single = nextrestart;
        modify->addstep_compute(next_restart_single);
      }
    }
    if (next_restart_double == ntimestep) {
      if (last_restart != ntimestep) {
        if (restart_toggle == 0) {
          restart->write(restart2a);
          restart_toggle = 1;
        } else {
          restart->write(restart2b);
          restart_toggle = 0;
        }
      }
      if (restart_every_double) next_restart_double += restart_every_double;
      else {
        modify->clearstep_compute();
        bigint nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_double));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad timestep");
        next_restart_double = nextrestart;
        modify->addstep_compute(next_restart_double);
      }
    }
    last_restart = ntimestep;
    next_restart = MIN(next_restart_single,next_restart_double);
  }

  // insure next_thermo forces output on last step of run
  // thermo may invoke computes so wrap with clear/add

  if (next_stats == ntimestep) {
    modify->clearstep_compute();
    if (last_stats != ntimestep) stats->compute(1);
    last_stats = ntimestep;
    if (var_stats) {
      next_stats = static_cast<bigint>
        (input->variable->compute_equal(ivar_stats));
      if (next_stats <= ntimestep)
        error->all(FLERR,"Stats every variable returned a bad timestep");
    } else if (stats_every) next_stats += stats_every;
    else next_stats = update->laststep;
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
   force restart file(s) to be written
------------------------------------------------------------------------- */

void Output::write_restart(bigint ntimestep)
{
  if (restart_flag_single) {
    char *file = new char[strlen(restart1) + 16];
    char *ptr = strchr(restart1,'*');
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s",restart1,ntimestep,ptr+1);
    *ptr = '*';
    restart->write(file);
    delete [] file;
  }

  if (restart_flag_double) {
    if (restart_toggle == 0) {
      restart->write(restart2a);
      restart_toggle = 1;
    } else {
      restart->write(restart2b);
      restart_toggle = 0;
    }
  }

  last_restart = ntimestep;
}

/* ----------------------------------------------------------------------
   timestep is being changed, called by update->reset_timestep()
   reset next timestep values for dumps, restart, thermo output
   reset to smallest value >= new timestep
   if next timestep set by variable evaluation,
     eval for ntimestep-1, so current ntimestep can be returned if needed
     no guarantee that variable can be evaluated for ntimestep-1
       if it depends on computes, but live with that rare case for now
------------------------------------------------------------------------- */

void Output::reset_timestep(bigint ntimestep)
{
  next_dump_any = MAXBIGINT;
  for (int idump = 0; idump < ndump; idump++) {
    if (every_dump[idump]) {
      next_dump[idump] = (ntimestep/every_dump[idump])*every_dump[idump];
      if (next_dump[idump] < ntimestep) next_dump[idump] += every_dump[idump];
    } else {
      modify->clearstep_compute();
      update->ntimestep--;
      bigint nextdump = static_cast<bigint>
        (input->variable->compute_equal(ivar_dump[idump]));
      if (nextdump < ntimestep)
        error->all(FLERR,"Dump every variable returned a bad timestep");
      update->ntimestep++;
      next_dump[idump] = nextdump;
      modify->addstep_compute(next_dump[idump]);
    }
    next_dump_any = MIN(next_dump_any,next_dump[idump]);
  }

  if (restart_flag_single) {
    if (restart_every_single) {
      next_restart_single =
        (ntimestep/restart_every_single)*restart_every_single;
      if (next_restart_single < ntimestep)
        next_restart_single += restart_every_single;
    } else {
      modify->clearstep_compute();
      update->ntimestep--;
      bigint nextrestart = static_cast<bigint>
        (input->variable->compute_equal(ivar_restart_single));
      if (nextrestart < ntimestep)
        error->all(FLERR,"Restart variable returned a bad timestep");
      update->ntimestep++;
      next_restart_single = nextrestart;
      modify->addstep_compute(next_restart_single);
    }
  } else next_restart_single = update->laststep + 1;

  if (restart_flag_double) {
    if (restart_every_double) {
      next_restart_double =
        (ntimestep/restart_every_double)*restart_every_double;
      if (next_restart_double < ntimestep)
        next_restart_double += restart_every_double;
    } else {
      modify->clearstep_compute();
      update->ntimestep--;
      bigint nextrestart = static_cast<bigint>
        (input->variable->compute_equal(ivar_restart_double));
      if (nextrestart < ntimestep)
        error->all(FLERR,"Restart variable returned a bad timestep");
      update->ntimestep++;
      next_restart_double = nextrestart;
      modify->addstep_compute(next_restart_double);
    }
  } else next_restart_double = update->laststep + 1;

  next_restart = MIN(next_restart_single,next_restart_double);

  if (var_stats) {
    modify->clearstep_compute();
    update->ntimestep--;
    next_stats = static_cast<bigint>
      (input->variable->compute_equal(ivar_stats));
    if (next_stats < ntimestep)
      error->all(FLERR,"Stats_modify every variable returned a bad timestep");
    update->ntimestep++;
    next_stats = MIN(next_stats,update->laststep);
    modify->addstep_compute(next_stats);
  } else if (stats_every) {
    next_stats = (ntimestep/stats_every)*stats_every;
    if (next_stats < ntimestep) next_stats += stats_every;
    next_stats = MIN(next_stats,update->laststep);
  } else next_stats = update->laststep;

  next = MIN(next_dump_any,next_restart);
  next = MIN(next,next_stats);
}

/* ----------------------------------------------------------------------
   add a Dump to list of Dumps
------------------------------------------------------------------------- */

void Output::add_dump(int narg, char **arg)
{
  if (narg < 5) error->all(FLERR,"Illegal dump command");

  // error checks

  for (int idump = 0; idump < ndump; idump++)
    if (strcmp(arg[0],dump[idump]->id) == 0)
      error->all(FLERR,"Reuse of dump ID");
  if (atoi(arg[3]) <= 0) error->all(FLERR,"Invalid dump frequency");

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
  else if (strcmp(arg[1],#key) == 0) dump[ndump] = new Class(sparta,narg,arg);
#include "style_dump.h"
#undef DUMP_CLASS

  else error->all(FLERR,"Unrecognized dump style");

  every_dump[ndump] = atoi(arg[3]);
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
   set stats output frequency from input script
------------------------------------------------------------------------- */

void Output::set_stats(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal stats command");

  if (strstr(arg[0],"v_") == arg[0]) {
    delete [] var_stats;
    int n = strlen(&arg[0][2]) + 1;
    var_stats = new char[n];
    strcpy(var_stats,&arg[0][2]);
  } else {
    stats_every = atoi(arg[0]);
    if (stats_every < 0) error->all(FLERR,"Illegal stats command");
  }
}

/* ----------------------------------------------------------------------
   new Stats style
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

  int every = 0;
  int varflag = 0;

  if (strstr(arg[0],"v_") == arg[0]) varflag = 1;
  else every = atoi(arg[0]);

  if (!varflag && every == 0) {
    if (narg != 1) error->all(FLERR,"Illegal restart command");

    restart_flag = restart_flag_single = restart_flag_double = 0;
    last_restart = -1;

    delete restart;
    restart = NULL;
    delete [] restart1;
    delete [] restart2a;
    delete [] restart2b;
    restart1 = restart2a = restart2b = NULL;
    delete [] var_restart_single;
    delete [] var_restart_double;
    var_restart_single = var_restart_double = NULL;

    return;
  }

  if (narg < 2) error->all(FLERR,"Illegal restart command");

  int nfile = 0;
  if (narg % 2 == 0) nfile = 1;
  else nfile = 2;

  if (nfile == 1) {
    restart_flag = restart_flag_single = 1;

    if (varflag) {
      delete [] var_restart_single;
      int n = strlen(&arg[0][2]) + 1;
      var_restart_single = new char[n];
      strcpy(var_restart_single,&arg[0][2]);
      restart_every_single = 0;
    } else restart_every_single = every;

    int n = strlen(arg[1]) + 3;
    restart1 = new char[n];
    strcpy(restart1,arg[1]);
    if (strchr(restart1,'*') == NULL) strcat(restart1,".*");
  }

  if (nfile == 2) {
    restart_flag = restart_flag_double = 1;

    if (varflag) {
      delete [] var_restart_double;
      int n = strlen(&arg[0][2]) + 1;
      var_restart_double = new char[n];
      strcpy(var_restart_double,&arg[0][2]);
      restart_every_double = 0;
    } else restart_every_double = every;

    restart_toggle = 0;
    int n = strlen(arg[1]) + 3;
    restart2a = new char[n];
    strcpy(restart2a,arg[1]);
    n = strlen(arg[2]) + 1;
    restart2b = new char[n];
    strcpy(restart2b,arg[2]);
  }

  // check for multiproc output and an MPI-IO filename
  // if 2 filenames, must be consistent

  int multiproc;
  if (strchr(arg[1],'%')) multiproc = comm->nprocs;
  else multiproc = 0;
  if (nfile == 2) {
    if (multiproc && !strchr(arg[2],'%'))
      error->all(FLERR,"Both restart files must use % or neither");
    if (!multiproc && strchr(arg[2],'%'))
      error->all(FLERR,"Both restart files must use % or neither");
  }

  // setup output style and process optional args

  delete restart;
  restart = new WriteRestart(sparta);
  int iarg = nfile+1;
  restart->multiproc_options(multiproc,narg-iarg,&arg[iarg]);
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

  MPI_Allreduce(&pbytes,&ave,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  double pave = scale * ave/comm->nprocs;
  MPI_Allreduce(&pbytes,&min,1,MPI_SPARTA_BIGINT,MPI_MIN,world);
  double pmin = scale * min;
  MPI_Allreduce(&pbytes,&max,1,MPI_SPARTA_BIGINT,MPI_MAX,world);
  double pmax = scale * max;

  MPI_Allreduce(&gbytes,&ave,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  double gave = scale * ave/comm->nprocs;
  MPI_Allreduce(&gbytes,&min,1,MPI_SPARTA_BIGINT,MPI_MIN,world);
  double gmin = scale * min;
  MPI_Allreduce(&gbytes,&max,1,MPI_SPARTA_BIGINT,MPI_MAX,world);
  double gmax = scale * max;

  MPI_Allreduce(&sbytes,&ave,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  double save = scale * ave/comm->nprocs;
  MPI_Allreduce(&sbytes,&min,1,MPI_SPARTA_BIGINT,MPI_MIN,world);
  double smin = scale * min;
  MPI_Allreduce(&sbytes,&max,1,MPI_SPARTA_BIGINT,MPI_MAX,world);
  double smax = scale * max;

  MPI_Allreduce(&bytes,&ave,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  double tave = scale * ave/comm->nprocs;
  MPI_Allreduce(&bytes,&min,1,MPI_SPARTA_BIGINT,MPI_MIN,world);
  double tmin = scale * min;
  MPI_Allreduce(&bytes,&max,1,MPI_SPARTA_BIGINT,MPI_MAX,world);
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
