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

#include "spapython.h"
#if defined(SPARTA_PYTHON)
#include "python_impl.h"
#else
#include "error.h"
#endif

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

Python::Python(SPARTA *sparta) : Pointers(sparta)
{
  // implementation of Python interface is only loaded on demand
  // and only if PYTHON package has been installed and compiled into binary
  impl = nullptr;
}

/* ---------------------------------------------------------------------- */

Python::~Python()
{
  delete impl;
}

/* ---------------------------------------------------------------------- */

void Python::init()
{
#if defined(SPARTA_PYTHON)
  if (!impl) impl = new PythonImpl(sparta);
#else
  error->all(FLERR, "Python support missing! Compile with PYTHON package installed!");
#endif
}

/* ---------------------------------------------------------------------- */
bool Python::is_enabled() const
{
#if defined(SPARTA_PYTHON)
  return true;
#else
  return false;
#endif
}

/* ---------------------------------------------------------------------- */

void Python::command(int narg, char **arg)
{
  init();
  impl->command(narg, arg);
}

/* ------------------------------------------------------------------ */

void Python::invoke_function(int ifunc, char *result, double *dvalue)
{
  init();
  impl->invoke_function(ifunc, result, dvalue);
}

/* ------------------------------------------------------------------ */

int Python::find(const char *name)
{
  init();
  return impl->find(name);
}

/* ------------------------------------------------------------------ */

int Python::function_match(const char *name, const char *varname, int numeric)
{
  init();
  return impl->function_match(name, varname, numeric);
}

/* ------------------------------------------------------------------ */

int Python::wrapper_match(const char *name, const char *varname,
                          int narg, int *argvars)
{
  init();
  return impl->wrapper_match(name, varname, narg, argvars);
}

/* ------------------------------------------------------------------ */

char *Python::long_string(int ifunc)
{
  init();
  return impl->long_string(ifunc);
}

/* ------------------------------------------------------------------ */

int Python::execute_string(char *cmd)
{
  init();
  return impl->execute_string(cmd);
}

/* ------------------------------------------------------------------ */

int Python::execute_file(char *fname)
{
  init();
  return impl->execute_file(fname);
}

/* ------------------------------------------------------------------ */

bool Python::has_minimum_version(int major, int minor)
{
  init();
  return impl->has_minimum_version(major, minor);
}

/* ------------------------------------------------------------------ */

void Python::finalize()
{
#if defined(SPARTA_PYTHON)
  PythonImpl::finalize();
#endif
}
