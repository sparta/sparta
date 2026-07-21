/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "error.h"
#include "spaexception.h"
#include "input.h"
#include "universe.h"
#include "output.h"
#include "memory.h"
#include "accelerator_kokkos.h"

using namespace SPARTA_NS;

// helper to format "message (file:line)" like the printed error text

static std::string truncpath(const char *path)
{
  if (path) {
    std::string full(path);
    std::size_t src = full.rfind("src/");
    if (src != std::string::npos) return full.substr(src);
    return full;
  }
  return "(unknown)";
}

static std::string fmt_error(const char *str, const char *file, int line)
{
  std::string msg = str ? str : "(unknown error)";
  msg += " (";
  msg += truncpath(file);
  msg += ":";
  msg += std::to_string(line);
  msg += ")";
  return msg;
}

// the input-script line being processed when the error occurred, so the
// message can echo it as "Last command: ..." (matches the LAMMPS behavior)

static std::string last_command(Input *input)
{
  if (input && input->line && strlen(input->line)) return input->line;
  return "(unknown)";
}


/* ---------------------------------------------------------------------- */

Error::Error(SPARTA *sparta) : Pointers(sparta) {}

/* ----------------------------------------------------------------------
   called by all procs in universe
   write error message to universe screen and logfile
   throw SpartaException, which all procs in the universe unwind to,
   either main() which terminates or the library interface which recovers
------------------------------------------------------------------------- */

void Error::universe_all(const char *file, int line, const char *str)
{
  MPI_Barrier(universe->uworld);

  if (universe->me == 0) {
    if (universe->uscreen) fprintf(universe->uscreen,
                                   "ERROR: %s (%s:%d)\n",str,file,line);
    if (universe->ulogfile) fprintf(universe->ulogfile,
                                    "ERROR: %s (%s:%d)\n",str,file,line);
    if (universe->uscreen) fflush(universe->uscreen);
    if (universe->ulogfile) fflush(universe->ulogfile);
  }

  std::string msg = fmt_error(str,file,line);
  set_last_error(msg.c_str(),ERROR_NORMAL);
  throw SpartaException(msg);
}

/* ----------------------------------------------------------------------
   called by one proc in universe
   throw SpartaAbortException; catch site must MPI_Abort if parallel,
   can recover if serial
------------------------------------------------------------------------- */

void Error::universe_one(const char *file, int line, const char *str)
{
  if (universe->uscreen) {
    fprintf(universe->uscreen,"ERROR on proc %d: %s (%s:%d)\n",
            universe->me,str,file,line);
    fflush(universe->uscreen);
  }

  std::string msg = fmt_error(str,file,line);
  set_last_error(msg.c_str(),ERROR_ABORT);
  throw SpartaAbortException(msg,universe->uworld);
}

/* ----------------------------------------------------------------------
   called by all procs in one world
   write error message to world screen and logfile
   throw SpartaException, which all procs in the world unwind to,
   either main() which terminates or the library interface which recovers
------------------------------------------------------------------------- */

void Error::all(const char *file, int line, const char *str)
{
  MPI_Barrier(world);

  int me;
  MPI_Comm_rank(world,&me);

  std::string lastcmd = last_command(input);

  if (me == 0) {
    if (screen) {
      fprintf(screen,"ERROR: %s (%s:%d)\n",str,file,line);
      fprintf(screen,"Last command: %s\n",lastcmd.c_str());
      fflush(screen);
    }
    if (logfile) {
      fprintf(logfile,"ERROR: %s (%s:%d)\n",str,file,line);
      fprintf(logfile,"Last command: %s\n",lastcmd.c_str());
      fflush(logfile);
    }
  }

  std::string msg = fmt_error(str,file,line);
  msg += "\nLast command: " + lastcmd;
  set_last_error(msg.c_str(),ERROR_NORMAL);
  throw SpartaException(msg);
}

/* ----------------------------------------------------------------------
   called by one proc in world
   write to world screen only if non-NULL on this proc
   always write to universe screen
   throw SpartaAbortException; catch site must MPI_Abort if parallel,
   can recover if serial
------------------------------------------------------------------------- */

void Error::one(const char *file, int line, const char *str)
{
  int me;
  MPI_Comm_rank(world,&me);

  std::string lastcmd = last_command(input);

  if (screen) {
    fprintf(screen,"ERROR on proc %d: %s (%s:%d)\n",
            me,str,file,line);
    fprintf(screen,"Last command: %s\n",lastcmd.c_str());
    fflush(screen);
  }
  if (universe->nworlds > 1) {
    fprintf(universe->uscreen,"ERROR on proc %d: %s (%s:%d)\n",
            universe->me,str,file,line);
    fprintf(universe->uscreen,"Last command: %s\n",lastcmd.c_str());
    fflush(universe->uscreen);
  }

  std::string msg = fmt_error(str,file,line);
  msg += "\nLast command: " + lastcmd;
  set_last_error(msg.c_str(),ERROR_ABORT);
  throw SpartaAbortException(msg,world);
}

/* ----------------------------------------------------------------------
   called by one proc in world
   only write to screen if non-NULL on this proc since could be file
------------------------------------------------------------------------- */

void Error::warning(const char *file, int line, const char *str, int logflag)
{
  if (screen) fprintf(screen,"WARNING: %s (%s:%d)\n",str,file,line);
  if (logflag && logfile) fprintf(logfile,"WARNING: %s (%s:%d)\n",
                                  str,file,line);
}

/* ----------------------------------------------------------------------
   called by one proc in world, typically proc 0
   write message to screen and logfile (if logflag is set)
------------------------------------------------------------------------- */

void Error::message(const char *file, int line, const char *str, int logflag)
{
  if (screen) fprintf(screen,"%s (%s:%d)\n",str,file,line);
  if (logflag && logfile) fprintf(logfile,"%s (%s:%d)\n",str,file,line);
}

/* ----------------------------------------------------------------------
   called by all procs in one world
   close all output, screen, and log files in world
   this terminates the process and is only invoked by the quit command
------------------------------------------------------------------------- */

void Error::done()
{
  MPI_Barrier(world);

  if (output) delete output;
  if (screen && screen != stdout) fclose(screen);
  if (logfile) fclose(logfile);

  if (sparta->kokkos) Kokkos::finalize();
  MPI_Finalize();
  exit(1);
}

/* ----------------------------------------------------------------------
   store the last error message and its type
   for retrieval via the library interface
------------------------------------------------------------------------- */

void Error::set_last_error(const char *msg, int type)
{
  if (msg) {
    last_error_message = msg;
    last_error_type = type;
  } else {
    last_error_message.clear();
    last_error_type = ERROR_NONE;
  }
}
