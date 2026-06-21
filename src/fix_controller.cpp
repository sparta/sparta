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
   File adapted from LAMMPS (https://www.lammps.org), October 2024
   Ported to SPARTA by: Stan Moore (SNL)
   Original Author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "fix_controller.h"
#include "compute.h"
#include "fix.h"
#include "error.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

using namespace SPARTA_NS;

enum { COMPUTE, FIX, VARIABLE };

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2

/* ---------------------------------------------------------------------- */

FixController::FixController(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg),
  pvID(nullptr), cvID(nullptr)
{
  if (narg != 10) error->all(FLERR,"Illegal fix controller command");

  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;

  nevery = atoi(arg[2]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix controller command");

  alpha = atof(arg[3]);
  kp = atof(arg[4]);
  ki = atof(arg[5]);
  kd = atof(arg[6]);

  // process variable arg = arg[7] : c_ID, c_ID[N], f_ID, f_ID[N], or v_ID

  if (strncmp(arg[7],"c_",2) == 0) pvwhich = COMPUTE;
  else if (strncmp(arg[7],"f_",2) == 0) pvwhich = FIX;
  else if (strncmp(arg[7],"v_",2) == 0) pvwhich = VARIABLE;
  else error->all(FLERR,"Illegal fix controller command");

  int n = strlen(arg[7]);
  char *suffix = new char[n];
  strcpy(suffix,&arg[7][2]);

  char *ptr = strchr(suffix,'[');
  if (ptr) {
    if (suffix[strlen(suffix)-1] != ']')
      error->all(FLERR,"Illegal fix controller command");
    pvindex = atoi(ptr+1);
    *ptr = '\0';
  } else pvindex = 0;

  n = strlen(suffix) + 1;
  pvID = new char[n];
  strcpy(pvID,suffix);
  delete [] suffix;

  // setpoint arg

  setpoint = atof(arg[8]);

  // control variable arg

  n = strlen(arg[9]) + 1;
  cvID = new char[n];
  strcpy(cvID,arg[9]);

  // error checks

  if (pvwhich == COMPUTE) {
    int icompute = modify->find_compute(pvID);
    if (icompute < 0)
      error->all(FLERR,"Compute ID for fix controller does not exist");
    Compute *c = modify->compute[icompute];
    int flag = 0;
    if (c->scalar_flag && pvindex == 0) flag = 1;
    else if (c->vector_flag && pvindex > 0) flag = 1;
    if (!flag)
      error->all(FLERR,"Fix controller compute does not calculate a "
                 "global scalar or vector");
    if (pvindex && pvindex > c->size_vector)
      error->all(FLERR,"Fix controller compute vector is accessed out-of-range");

  } else if (pvwhich == FIX) {
    int ifix = modify->find_fix(pvID);
    if (ifix < 0)
      error->all(FLERR,"Fix ID for fix controller does not exist");
    Fix *f = modify->fix[ifix];
    int flag = 0;
    if (f->scalar_flag && pvindex == 0) flag = 1;
    else if (f->vector_flag && pvindex > 0) flag = 1;
    if (!flag)
      error->all(FLERR,"Fix controller fix does not calculate a "
                 "global scalar or vector");
    if (pvindex && pvindex > f->size_vector)
      error->all(FLERR,"Fix controller fix vector is accessed out-of-range");

  } else if (pvwhich == VARIABLE) {
    int ivariable = input->variable->find(pvID);
    if (ivariable < 0)
      error->all(FLERR,"Variable name for fix controller does not exist");
    if (input->variable->equal_style(ivariable) == 0)
      error->all(FLERR,"Fix controller variable is not equal-style variable");
  }

  int ivariable = input->variable->find(cvID);
  if (ivariable < 0)
    error->all(FLERR,"Variable name for fix controller does not exist");
  if (input->variable->internal_style(ivariable) == 0)
    error->all(FLERR,"Fix controller variable is not internal-style variable");
  control = input->variable->compute_equal(ivariable);

  firsttime = 1;
}

/* ---------------------------------------------------------------------- */

FixController::~FixController()
{
  delete [] pvID;
  delete [] cvID;
}

/* ---------------------------------------------------------------------- */

int FixController::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixController::init()
{
  if (pvwhich == COMPUTE) {
    int icompute = modify->find_compute(pvID);
    if (icompute < 0)
      error->all(FLERR,"Compute ID for fix controller does not exist");
    pcompute = modify->compute[icompute];

  } else if (pvwhich == FIX) {
    int ifix = modify->find_fix(pvID);
    if (ifix < 0)
      error->all(FLERR,"Fix ID for fix controller does not exist");
    pfix = modify->fix[ifix];

  } else if (pvwhich == VARIABLE) {
    pvar = input->variable->find(pvID);
    if (pvar < 0)
      error->all(FLERR,"Variable name for fix controller does not exist");
  }

  cvar = input->variable->find(cvID);
  if (cvar < 0)
    error->all(FLERR,"Variable name for fix controller does not exist");

  // set sampling time

  tau = nevery * update->dt;
}

/* ---------------------------------------------------------------------- */

void FixController::end_of_step()
{
  // current value of pv = invocation of compute,fix,variable
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  double current = 0.0;

  if (pvwhich == COMPUTE) {

    // invoke compute if not previously invoked

    if (pvindex == 0) {
      if (!(pcompute->invoked_flag & INVOKED_SCALAR)) {
        pcompute->compute_scalar();
        pcompute->invoked_flag |= INVOKED_SCALAR;
      }
      current = pcompute->scalar;
    } else {
      if (!(pcompute->invoked_flag & INVOKED_VECTOR)) {
        pcompute->compute_vector();
        pcompute->invoked_flag |= INVOKED_VECTOR;
      }
      current = pcompute->vector[pvindex-1];
    }

  // access fix field, guaranteed to be ready

  } else if (pvwhich == FIX) {
    if (pvindex == 0) current = pfix->compute_scalar();
    else current = pfix->compute_vector(pvindex-1);

  // evaluate equal-style variable

  } else if (pvwhich == VARIABLE) {
    current = input->variable->compute_equal(pvar);
  }

  modify->addstep_compute(update->ntimestep + nevery);

  // new control var = f(old value, current process var, setpoint)
  // cv = cvold -kp*err -ki*sumerr -kd*deltaerr
  // note: this deviates from standard notation, which is
  // cv = kp*err +ki*sumerr +kd*deltaerr
  // the difference is in the sign and the time integral

  err = current - setpoint;

  if (firsttime) {
    firsttime = 0;
    deltaerr = sumerr = 0.0;
  } else {
    deltaerr = err - olderr;
  }
  sumerr += err;

  // 3 terms of PID equation

  control += -kp * alpha * tau * err;
  control += -ki * alpha * tau * tau * sumerr;
  control += -kd * alpha * deltaerr;
  olderr = err;

  // reset control variable

  input->variable->internal_set(cvar,control);
}

/* ----------------------------------------------------------------------
   return 3 terms of PID controller at last invocation of end_of_step()
------------------------------------------------------------------------- */

double FixController::compute_vector(int n)
{
  if (n == 0) return (-kp * alpha * tau * err);
  else if (n == 1) return (-ki * alpha * tau * tau * sumerr);
  else return (-kd * alpha * deltaerr);
}
