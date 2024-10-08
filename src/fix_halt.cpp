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

/* ----------------------------------------------------------------------
   File copied and adapted from LAMMPS (lammps.org), October 2024
   Author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "fix_halt.h"

#include "comm.h"
#include "error.h"
#include "input.h"
#include "modify.h"
#include "timer.h"
#include "update.h"
#include "utils.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace SPARTA_NS;

enum { TLIMIT, VARIABLE };
enum { LT, LE, GT, GE, EQ, NEQ, XOR };
enum { HARD, SOFT, CONTINUE };
enum { NOMSG = 0, YESMSG = 1 };

/* ---------------------------------------------------------------------- */

FixHalt::FixHalt(SPARTA *sparta, int narg, char **arg) :
    Fix(sparta, narg, arg), idvar(nullptr)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix halt", error);
  nevery = utils::inumeric(FLERR, arg[2], false, sparta);
  if (nevery <= 0) error->all(FLERR, "Illegal fix halt command: nevery must be > 0");

  // comparison args

  idvar = nullptr;
  int iarg = 3;

  if (strcmp(arg[iarg], "tlimit") == 0) {
    attribute = TLIMIT;
  } else {
    if (!utils::strmatch(arg[iarg],"^v_")) {
      char msg[128];
      sprintf(msg, "Invalid fix halt attribute %s", arg[iarg]);
      error->all(FLERR, msg);
    }

    char* str = arg[iarg];
    attribute = VARIABLE;
    int n = strlen(&str[2]) + 1;
    idvar = new char[n];
    strcpy(idvar,&str[2]);
    ivar = input->variable->find(idvar);

    if (ivar < 0) error->all(FLERR, "Could not find fix halt variable name");
    if (input->variable->equal_style(ivar) == 0)
      error->all(FLERR, "Fix halt variable is not equal-style variable");
  }

  // clang-format off
  ++iarg;
  if (strcmp(arg[iarg],"<") == 0) operation = LT;
  else if (strcmp(arg[iarg],"<=") == 0) operation = LE;
  else if (strcmp(arg[iarg],">") == 0) operation = GT;
  else if (strcmp(arg[iarg],">=") == 0) operation = GE;
  else if (strcmp(arg[iarg],"==") == 0) operation = EQ;
  else if (strcmp(arg[iarg],"!=") == 0) operation = NEQ;
  else if (strcmp(arg[iarg],"|^") == 0) operation = XOR;
  else error->all(FLERR,"Invalid fix halt operator");

  ++iarg;
  value = utils::numeric(FLERR, arg[iarg], false, sparta);

  // parse optional args

  eflag = SOFT;
  msgflag = YESMSG;
  ++iarg;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "error") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix halt error", error);
      if (strcmp(arg[iarg + 1], "hard") == 0) eflag = HARD;
      else if (strcmp(arg[iarg + 1], "soft") == 0) eflag = SOFT;
      else if (strcmp(arg[iarg + 1], "continue") == 0) eflag = CONTINUE;
      else {
        char msg[128];
        sprintf(msg, "Unknown fix halt error condition %s", arg[iarg]);
        error->all(FLERR, msg);
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "message") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix halt message", error);
      msgflag = utils::logical(FLERR, arg[iarg + 1], false, sparta);
      iarg += 2;
    } else {
      char msg[128];
      sprintf(msg, "Unknown fix halt keyword %s", arg[iarg]);
      error->all(FLERR, msg);
    }
  }
  // clang-format on

  // add nfirst to all computes that store invocation times
  // since don't know a priori which are invoked via variables by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  if (attribute == VARIABLE) {
    const bigint nfirst = (update->ntimestep / nevery) * nevery + nevery;
    modify->addstep_compute_all(nfirst);
  }
}

/* ---------------------------------------------------------------------- */

FixHalt::~FixHalt()
{
  delete[] idvar;
}

/* ---------------------------------------------------------------------- */

int FixHalt::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= POST_RUN;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHalt::init()
{
  // set ivar from current variable list

  if (attribute == VARIABLE) {
    ivar = input->variable->find(idvar);
    char msg[128];
    if (ivar < 0) {
      sprintf(msg, "Could not find fix halt variable %s", idvar);
      error->all(FLERR, msg);
    }
    if (input->variable->equal_style(ivar) == 0) {
      sprintf(msg, "Fix halt variable %s is not equal-style variable", idvar);
      error->all(FLERR, msg);
    }
  }

  // settings used by TLIMIT

  nextstep = (update->ntimestep / nevery) * nevery + nevery;
  thisstep = -1;
  tratio = 0.5;
}

/* ---------------------------------------------------------------------- */

void FixHalt::end_of_step()
{
  // variable evaluation may invoke computes so wrap with clear/add

  double attvalue;

  if (attribute == TLIMIT) {
    if (update->ntimestep != nextstep) return;
    attvalue = tlimit();
  } else {
    modify->clearstep_compute();
    attvalue = input->variable->compute_equal(ivar);
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // ensure that the attribute is *exactly* the same on all ranks

  MPI_Bcast(&attvalue, 1, MPI_DOUBLE, 0, world);

  // check if halt is triggered, else just return

  if (operation == LT) {
    if (attvalue >= value) return;
  } else if (operation == LE) {
    if (attvalue > value) return;
  } else if (operation == GT) {
    if (attvalue <= value) return;
  } else if (operation == GE) {
    if (attvalue < value) return;
  } else if (operation == EQ) {
    if (attvalue != value) return;
  } else if (operation == NEQ) {
    if (attvalue == value) return;
  } else if (operation == XOR) {
    if ((attvalue == 0.0 && value == 0.0) || (attvalue != 0.0 && value != 0.0)) return;
  }

  // hard halt -> exit SPARTA
  // soft/continue halt -> trigger timer to break from run loop
  // print message with ID of fix halt in case multiple instances

  char message[128];
  sprintf(message, "Fix halt condition for fix-id %s met on step %ld with value %g",
                                    id, update->ntimestep, attvalue);
  if (eflag == HARD) {
    error->all(FLERR, message);
  } else if ((eflag == SOFT) || (eflag == CONTINUE)) {
    if ((comm->me == 0) && (msgflag == YESMSG)) error->message(FLERR, message);
    timer->force_timeout();
  }
}

/* ----------------------------------------------------------------------
   reset expired timer setting to original value, if requested
------------------------------------------------------------------------- */

void FixHalt::post_run()
{
  // continue halt -> subsequent runs are allowed

  if (eflag == CONTINUE) timer->reset_timeout();
}

/* ----------------------------------------------------------------------
   compute synced elapsed time
   reset nextstep = estimate of timestep when run will end
   first project to 1/2 the run time, thereafter to end of run
------------------------------------------------------------------------- */

double FixHalt::tlimit()
{
  double cpu = timer->elapsed(TIME_LOOP);
  MPI_Bcast(&cpu, 1, MPI_DOUBLE, 0, world);

  if (cpu < value) {
    bigint elapsed = update->ntimestep - update->firststep;
    bigint final = update->firststep + static_cast<bigint>(tratio * value / cpu * elapsed);
    nextstep = (final / nevery) * nevery + nevery;
    if (nextstep == update->ntimestep) nextstep += nevery;
    tratio = 1.0;
  }

  return cpu;
}
