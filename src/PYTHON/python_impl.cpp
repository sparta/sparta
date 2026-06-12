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
   Contributing authors: Richard Berger (LANL) and Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "python_impl.h"

#include "comm.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "python_compat.h"
#include "python_utils.h"
#include "utils.h"
#include "variable.h"

#include <Python.h>    // IWYU pragma: export
#include <cstring>

using namespace SPARTA_NS;

enum { NONE, INT, DOUBLE, STRING, PTR };
enum { VALUE, VARIABLE, INTERNALVAR };

/* ---------------------------------------------------------------------- */

PythonImpl::PythonImpl(SPARTA *sparta) : Pointers(sparta)
{
  // pfuncs stores interface info for each Python function

  nfunc = 0;
  pfuncs = nullptr;

#if PY_MAJOR_VERSION >= 3 && !defined(Py_LIMITED_API)
  // check for PYTHONUNBUFFERED environment variable
  const char *PYTHONUNBUFFERED = getenv("PYTHONUNBUFFERED");
  // Force the stdout and stderr streams to be unbuffered.
  bool unbuffered = PYTHONUNBUFFERED != nullptr && strcmp(PYTHONUNBUFFERED, "1") == 0;

#if (PY_VERSION_HEX >= 0x030800f0)
  PyConfig config;
  PyConfig_InitPythonConfig(&config);
  config.buffered_stdio = !unbuffered;
#else
  // Python Global configuration variable
  Py_UnbufferedStdioFlag = unbuffered;
#endif
#endif

#if PY_VERSION_HEX >= 0x030800f0 && !defined(Py_LIMITED_API)
  Py_InitializeFromConfig(&config);
  PyConfig_Clear(&config);
#else
  Py_Initialize();
#endif

  // only needed for Python 2.x and Python 3 < 3.7
  // With Python 3.7 this function is now called by Py_Initialize()
  // Deprecated since version 3.9, will be removed in version 3.11
#if PY_VERSION_HEX < 0x030700f0
  if (!PyEval_ThreadsInitialized()) PyEval_InitThreads();
#endif

  PyUtils::GIL lock;

  PyObject *pModule = PyImport_AddModule("__main__");
  if (!pModule) error->all(FLERR, "Could not initialize embedded Python");

  pyMain = (void *) pModule;
}

/* ---------------------------------------------------------------------- */

PythonImpl::~PythonImpl()
{
  if (pyMain) {
    // clean up
    PyUtils::GIL lock;

    for (int i = 0; i < nfunc; i++) {
      delete[] pfuncs[i].name;
      deallocate(i);
      Py_CLEAR(pfuncs[i].pFunc);
    }
  }

  memory->sfree(pfuncs);
}

/* ---------------------------------------------------------------------- */

void PythonImpl::command(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal python command");

  // if invoke is only keyword, invoke the previously defined function

  if (strcmp(arg[1], "invoke") == 0) {
    int ifunc = find(arg[0]);
    if (ifunc < 0) {
      char msg[128];
      sprintf(msg, "Python invoke of unknown function: %s", arg[0]);
      error->all(FLERR, msg);
    }

    char *str = nullptr;
    if (pfuncs[ifunc].noutput) {
      str = input->variable->python_style(pfuncs[ifunc].ovarname, pfuncs[ifunc].name);
      if (!str) {
        char msg[128];
        sprintf(msg,
                   "Python variable %s does not match variable %s "
                   "registered with Python function %s",
                   arg[0], pfuncs[ifunc].ovarname, pfuncs[ifunc].name);
        error->all(FLERR, msg);
      }
    }

    bool logreturn = false;
    if (narg == 3 && strcmp(arg[2], "logreturn") == 0) logreturn = true;

    invoke_function(ifunc, str, NULL);

    if (logreturn && str && (comm->me == 0)) {
      if (screen) fprintf(screen,"Invoked python function %s returned %s\n",arg[0],str);
      if (logfile) fprintf(logfile,"Invoked python function %s returned %s\n",arg[0],str);
    }

    return;
  }

  // if source is only keyword, execute the python code in file

  if ((narg > 1) && (strcmp(arg[0], "source") == 0)) {
    int err = -1;

    if ((narg > 2) && (strcmp(arg[1], "here") == 0)) {
      err = execute_string(arg[2]);
    } else {
      int file_is_readable = 0;
      FILE *fp = fopen(arg[1], "r");
      if (fp) {
        fclose(fp);
        file_is_readable = 1;
      }
      if (file_is_readable)
        err = execute_file(arg[1]);
      else {
        char msg[128];
        sprintf(msg, "Could not open python source file %s for processing", arg[1]);
        error->all(FLERR, msg);
      }
    }
    if (err) error->all(FLERR, "Failure in python source command");

    return;
  }

  // parse optional args, invoke is not allowed in this mode

  int ninput = 0;
  int noutput = 0;
  char **istr = nullptr;
  char *ostr = nullptr;
  char *format = nullptr;
  int length_longstr = 0;
  char *pyfile = nullptr;
  char *herestr = nullptr;
  int existflag = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "input") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "python input", error);
      ninput = utils::inumeric(FLERR, arg[iarg + 1], false, sparta);
      if (ninput < 0) {
        char msg[128];
        sprintf(msg, "Invalid number of python input arguments: %i", ninput);
        error->all(FLERR, msg);
      }
      iarg += 2;
      delete[] istr;
      istr = new char *[ninput];
      if (iarg+ninput > narg) error->all(FLERR, "Invalid python input command");
      for (int i = 0; i < ninput; i++) istr[i] = arg[iarg + i];
      iarg += ninput;
    } else if (strcmp(arg[iarg], "return") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Invalid python return command");
      noutput = 1;
      ostr = arg[iarg + 1];
      iarg += 2;
    } else if (strcmp(arg[iarg], "format") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Invalid python format command");
      format = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "length") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "python length", error);
      length_longstr = utils::inumeric(FLERR, arg[iarg + 1], false, sparta);
      if (length_longstr <= 0) error->all(FLERR, "Invalid python return value length");
      iarg += 2;
    } else if (strcmp(arg[iarg], "file") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Invalid python file command");
      delete[] pyfile;
      pyfile = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "here") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Invalid python here command");
      herestr = arg[iarg + 1];
      iarg += 2;
    } else if (strcmp(arg[iarg], "exists") == 0) {
      existflag = 1;
      iarg++;
    } else {
      char msg[128];
      sprintf(msg, "Unknown python command keyword: %s", arg[iarg]);
      error->all(FLERR, msg);
    }
  }

  if (pyfile && herestr)
    error->all(FLERR, "Must not use python 'file' and 'here' keywords at the same time");
  if (pyfile && existflag)
    error->all(FLERR, "Must not use python 'file' and 'exists' keywords at the same time");
  if (herestr && existflag)
    error->all(FLERR, "Must not use python 'here' and 'exists' keywords at the same time");

  // create or overwrite entry in pfuncs vector with name = arg[0]

  int ifunc = create_entry(arg[0], ninput, noutput, length_longstr, istr, ostr, format);

  PyUtils::GIL lock;

  // send Python code to Python interpreter
  // file: read the file via PyRun_SimpleFile()
  // here: process the here string directly
  // exist: do nothing, assume code has already been run

  if (pyfile) {
    FILE *fp = fopen(pyfile, "r");

    if (fp == nullptr) {
      PyUtils::Print_Errors();
      char msg[128];
      sprintf(msg, "Could not open Python file: %s", pyfile);
      error->all(FLERR, msg);
    }

    int err = PyRun_SimpleFile(fp, pyfile);
    if (err) {
      PyUtils::Print_Errors();
      char msg[128];
      sprintf(msg, "Could not open Python file: %s", pyfile);
      error->all(FLERR, msg);
    }
    fclose(fp);

  } else if (herestr) {
    int err = PyRun_SimpleString(herestr);
    if (err) {
      PyUtils::Print_Errors();
      char msg[128];
      sprintf(msg, "Could not process Python string: %s", herestr);
      error->all(FLERR, msg);
    }
  }

  // pFunc = function object for requested function

  auto pModule = (PyObject *) pyMain;
  PyObject *pFunc = PyObject_GetAttrString(pModule, pfuncs[ifunc].name);

  if (!pFunc) {
    PyUtils::Print_Errors();
    char msg[128];
    sprintf(msg, "Could not find Python function %s", pfuncs[ifunc].name);
    error->all(FLERR, msg);
  }

  if (!PyCallable_Check(pFunc)) {
    PyUtils::Print_Errors();
    char msg[128];
    sprintf(msg, "Python function %s is not callable", pfuncs[ifunc].name);
    error->all(FLERR, msg);
  }

  pfuncs[ifunc].pFunc = (void *) pFunc;

  // clean-up input storage

  delete[] istr;
  delete[] format;
  delete[] pyfile;
}

/* ------------------------------------------------------------------ */

void PythonImpl::invoke_function(int ifunc, char *result, double *dvalue)
{
  PyUtils::GIL lock;
  PyObject *pValue;
  char *str;

  auto pFunc = (PyObject *) pfuncs[ifunc].pFunc;

  // create Python tuple of input arguments

  int ninput = pfuncs[ifunc].ninput;
  PyObject *pArgs = PyTuple_New(ninput);

  if (!pArgs) {
    char msg[128];
    sprintf(msg, "Could not prepare arguments for Python function %s", pfuncs[ifunc].name);
    error->all(FLERR, msg);
  }

  for (int i = 0; i < ninput; i++) {
    int itype = pfuncs[ifunc].itype[i];
    if (itype == INT) {
      if (pfuncs[ifunc].ivarflag[i] == VARIABLE) {
        str = input->variable->retrieve(pfuncs[ifunc].svalue[i]);
        if (!str) {
          char msg[128];
          sprintf(msg, "Could not evaluate Python function %s input variable: %s",
                     pfuncs[ifunc].name, pfuncs[ifunc].svalue[i]);
          error->all(FLERR, msg);
        }
        pValue = PY_INT_FROM_LONG(PY_LONG_FROM_STRING(str));
      } else if (pfuncs[ifunc].ivarflag[i] == INTERNALVAR) {
        double value = input->variable->compute_equal(pfuncs[ifunc].internal_var[i]);
        pValue = PyLong_FromDouble(value);
      } else {
        pValue = PY_INT_FROM_LONG(pfuncs[ifunc].ivalue[i]);
      }
    } else if (itype == DOUBLE) {
      if (pfuncs[ifunc].ivarflag[i] == VARIABLE) {
        str = input->variable->retrieve(pfuncs[ifunc].svalue[i]);
        if (!str) {
          char msg[128];
          sprintf(msg, "Could not evaluate Python function %s input variable: %s",
                     pfuncs[ifunc].name, pfuncs[ifunc].svalue[i]);
          error->all(FLERR, msg);
        }
        pValue = PyFloat_FromDouble(std::stod(str));
      } else if (pfuncs[ifunc].ivarflag[i] == INTERNALVAR) {
        double value = input->variable->compute_equal(pfuncs[ifunc].internal_var[i]);
        pValue = PyFloat_FromDouble(value);
      } else {
        pValue = PyFloat_FromDouble(pfuncs[ifunc].dvalue[i]);
      }
    } else if (itype == STRING) {
      if (pfuncs[ifunc].ivarflag[i] == VARIABLE) {
        str = input->variable->retrieve(pfuncs[ifunc].svalue[i]);
        if (!str) {
          char msg[128];
          sprintf(msg, "Could not evaluate Python function %s input variable: %s",
                     pfuncs[ifunc].name, pfuncs[ifunc].svalue[i]);
          error->all(FLERR, msg);
        }
        pValue = PyUnicode_FromString(str);
      } else {
        pValue = PyUnicode_FromString(pfuncs[ifunc].svalue[i]);
      }
    } else if (itype == PTR) {
      pValue = PyCapsule_New((void *)sparta, nullptr, nullptr);
    } else {
      char msg[128];
      sprintf(msg, "Unsupported variable type: %i", itype);
      error->all(FLERR, msg);
    }
    PyTuple_SetItem(pArgs, i, pValue);
  }

  // call the Python function
  // error check with one() since only some procs may fail

  pValue = PyObject_CallObject(pFunc, pArgs);
  Py_CLEAR(pArgs);

  if (!pValue) {
    PyUtils::Print_Errors();
    char msg[128];
    sprintf(msg, "Python evaluation of function %s failed", pfuncs[ifunc].name);
    error->all(FLERR, msg);
  }

  // function returned a value
  // if result is non-NULL, assign to result string stored by python-style variable
  //   or if value is string and user specified a length, assign it to longstr
  // if dvalue is non-NULL, assign numeric value directly to dvalue

  char python2str[64];

  if (pfuncs[ifunc].noutput) {
    int otype = pfuncs[ifunc].otype;
    if (otype == INT) {
      if (dvalue) *dvalue = (double) PY_INT_AS_LONG(pValue);
      else {
        char value[128];
        sprintf(value, "%ld", PY_INT_AS_LONG(pValue));
        strncpy(result, value, Variable::VALUELENGTH - 1);
      }
    } else if (otype == DOUBLE) {
      if (dvalue) *dvalue = PyFloat_AsDouble(pValue);
      else {
        char value[128];
        sprintf(value, "%.15g", PyFloat_AsDouble(pValue));
        strncpy(result, value, Variable::VALUELENGTH - 1);
      }
    } else if (otype == STRING) {
      const char *pystr = PyUnicode_AsUTF8(pValue);
      if (pfuncs[ifunc].longstr)
        strncpy(pfuncs[ifunc].longstr, pystr, pfuncs[ifunc].length_longstr);
      else
        strncpy(result, pystr, Variable::VALUELENGTH - 1);
    }
  }
  Py_CLEAR(pValue);
}

/* ------------------------------------------------------------------ */

int PythonImpl::find(const char *name)
{
  for (int i = 0; i < nfunc; i++)
    if (strcmp(name, pfuncs[i].name) == 0) return i;
  return -1;
}

/* ---------------------------------------------------------------------
   called by Variable class when a python-style variable is evaluated
     this will call invoke_function() in this class
   either via Variable::retrieve() or Variable::equalstyle
     retrieve calls with numeric = 0, equalstyle with numeric = 1
   ensure name matches a Python function
   ensure the Python function produces an output
   ensure the Python function outputs to the matching python-style variable
   ensure a string is returned only if retrieve() is the caller
--------------------------------------------------------------------- */

int PythonImpl::function_match(const char *name, const char *varname, int numeric)
{
  int ifunc = find(name);

  if (ifunc < 0)
    error->all(FLERR,"Python function specified by variable not found");
  if (pfuncs[ifunc].noutput == 0)
    error->all(FLERR,"Python function for variable does not return a value");
  if (strcmp(pfuncs[ifunc].ovarname, varname) != 0)
    error->all(FLERR,"Python variable does not match variable name "
               "registered with Python function");
  if (numeric && pfuncs[ifunc].otype == STRING)
    error->all(FLERR,"Python function for variable returns a string");

  return ifunc;
}

/* ---------------------------------------------------------------------
   called by Variable class when evaluating a Python wrapper function
     which will call invoke_function()
   either via equal-style or atom-style variable formula
     the latter calls invoke_function() once per atom
   same error checks as function_match() plus 2 new ones
   ensure match of number of Python function args mapped to internal-style variables
   ensure each internal-style variable still exists
     must check now in case user input script deleted variables between runs
       which could invalidate indices set in create_event()
     other classes avoid this issue by setting variable indices in their init() method
--------------------------------------------------------------------- */

int PythonImpl::wrapper_match(const char *name, const char *varname, int narg, int *argvars)
{
  int ifunc = function_match(name,varname,1);

  int internal_count = 0;
  for (int i = 0; i < pfuncs[ifunc].ninput; i++)
    if (pfuncs[ifunc].ivarflag[i] == INTERNALVAR) internal_count++;
  if (internal_count != narg)
    error->all(FLERR, "Mismatch in Python function and variable "
               "use of internal variable args");

  // set argvars of internal-style variables for use by Variable class
  //   in Python wrapper functions
  // also set internal_var for use by invoke_function()
  //   so that invoke_function() is as fast as possible for args which are internal-style vars

  int j = 0;
  for (int i = 0; i < pfuncs[ifunc].ninput; i++)
    if (pfuncs[ifunc].ivarflag[i] == INTERNALVAR) {
      int ivar = input->variable->find(pfuncs[ifunc].svalue[i]);
      if (ivar < 0)
        error->all(FLERR,"Python function cannot find an internal variable");
      pfuncs[ifunc].internal_var[i] = ivar;
      argvars[j++] = ivar;
    }

  return ifunc;
}

/* ------------------------------------------------------------------ */

char *PythonImpl::long_string(int ifunc)
{
  return pfuncs[ifunc].longstr;
}

/* ------------------------------------------------------------------ */

int PythonImpl::create_entry(char *name, int ninput, int noutput,
                             int length_longstr, char **istr,
                             char *ostr, char *format)
{
  // ifunc = index to entry by name in pfuncs vector, can be old or new
  // free old vectors if overwriting old pfunc

  int ifunc = find(name);

  if (ifunc < 0) {
    ifunc = nfunc;
    nfunc++;
    pfuncs = (PyFunc *) memory->srealloc(pfuncs, nfunc * sizeof(struct PyFunc), "python:pfuncs");
    pfuncs[ifunc].name = utils::strdup(name);
  } else
    deallocate(ifunc);

  pfuncs[ifunc].ninput = ninput;
  pfuncs[ifunc].noutput = noutput;

  if (!format && ninput + noutput)
    error->all(FLERR, "Missing python format keyword");
  else if (format && ((int) strlen(format) != ninput + noutput)) {
    char msg[128];
    sprintf(msg, "Input/output arguments (%i) and format characters (%zu) are inconsistent",
               (ninput + noutput), strlen(format));
    error->all(FLERR, msg);
  }

  // process inputs as values or variables

  pfuncs[ifunc].itype = new int[ninput];
  pfuncs[ifunc].ivarflag = new int[ninput];
  pfuncs[ifunc].ivalue = new int[ninput];
  pfuncs[ifunc].dvalue = new double[ninput];
  pfuncs[ifunc].svalue = new char *[ninput];
  pfuncs[ifunc].internal_var = new int[ninput];

  for (int i = 0; i < ninput; i++) {
    pfuncs[ifunc].svalue[i] = nullptr;
    char type = format[i];
    if (type == 'i') {
      pfuncs[ifunc].itype[i] = INT;
      if (utils::strmatch(istr[i], "^v_")) {
        pfuncs[ifunc].ivarflag[i] = VARIABLE;
        pfuncs[ifunc].svalue[i] = utils::strdup(istr[i] + 2);
      } else if (utils::strmatch(istr[i], "^iv_")) {
        pfuncs[ifunc].ivarflag[i] = INTERNALVAR;
        pfuncs[ifunc].svalue[i] = utils::strdup(istr[i] + 3);
        char *vname = pfuncs[ifunc].svalue[i];
        int ivar = input->variable->find(vname);
        if (ivar < 0) {   // create internal variable if does not exist
          input->variable->internal_create(vname,0.0);
          ivar = input->variable->find(vname);
        }
        if (!input->variable->internal_style(ivar)) {
          char str[128];
          sprintf(str,"Variable %s for python command is invalid style",vname);
          error->all(FLERR, str);
        }
      } else {
        pfuncs[ifunc].ivarflag[i] = VALUE;
        pfuncs[ifunc].ivalue[i] = utils::inumeric(FLERR, istr[i], false, sparta);
      }
    } else if (type == 'f') {
      pfuncs[ifunc].itype[i] = DOUBLE;
      if (utils::strmatch(istr[i], "^v_")) {
        pfuncs[ifunc].ivarflag[i] = VARIABLE;
        pfuncs[ifunc].svalue[i] = utils::strdup(istr[i] + 2);
      } else if (utils::strmatch(istr[i], "^iv_")) {
        pfuncs[ifunc].ivarflag[i] = INTERNALVAR;
        pfuncs[ifunc].svalue[i] = utils::strdup(istr[i] + 3);
        char *vname = pfuncs[ifunc].svalue[i];
        int ivar = input->variable->find(vname);
        if (ivar < 0) {   // create internal variable if does not exist
          input->variable->internal_create(vname,0.0);
          ivar = input->variable->find(vname);
        }
        if (!input->variable->internal_style(ivar)) {
          char str[128];
          sprintf(str,"Variable %s for python command is invalid style",vname);
          error->all(FLERR, str);
        }
      } else {
        pfuncs[ifunc].ivarflag[i] = VALUE;
        pfuncs[ifunc].dvalue[i] = utils::numeric(FLERR, istr[i], false, sparta);
      }
    } else if (type == 's') {
      pfuncs[ifunc].itype[i] = STRING;
      if (utils::strmatch(istr[i], "^v_")) {
        pfuncs[ifunc].ivarflag[i] = VARIABLE;
        pfuncs[ifunc].svalue[i] = utils::strdup(istr[i] + 2);
      } else if (utils::strmatch(istr[i], "^iv_")) {
        char str[128];
        sprintf(str,"Input argument %s cannot be internal variable with string format",istr[i]);
        error->all(FLERR, str);
      } else {
        pfuncs[ifunc].ivarflag[i] = 0;
        pfuncs[ifunc].svalue[i] = utils::strdup(istr[i]);
      }
    } else if (type == 'p') {
      pfuncs[ifunc].ivarflag[i] = VALUE;
      pfuncs[ifunc].itype[i] = PTR;
      if (strcmp(istr[i], "SELF") != 0) error->all(FLERR, "Invalid python command");

    } else {
      char msg[128];
      sprintf(msg, "Invalid python format character: %i", type);
      error->all(FLERR, msg);
    }
  }

  // process output as value or variable

  pfuncs[ifunc].ovarname = nullptr;
  pfuncs[ifunc].longstr = nullptr;
  if (!noutput) return ifunc;

  char type = format[ninput];
  if (type == 'i')
    pfuncs[ifunc].otype = INT;
  else if (type == 'f')
    pfuncs[ifunc].otype = DOUBLE;
  else if (type == 's')
    pfuncs[ifunc].otype = STRING;
  else {
    char msg[128];
    sprintf(msg, "Invalid python return format character: %i", type);
    error->all(FLERR, msg);
  }

  if (length_longstr) {
    if (pfuncs[ifunc].otype != STRING)
      error->all(FLERR, "Python command length keyword cannot be used unless output is a string");
    pfuncs[ifunc].length_longstr = length_longstr;
    pfuncs[ifunc].longstr = new char[length_longstr + 1];
    pfuncs[ifunc].longstr[length_longstr] = '\0';
  }

  if (strstr(ostr, "v_") != ostr) error->all(FLERR, "Invalid python command");
  pfuncs[ifunc].ovarname = utils::strdup(ostr + 2);

  return ifunc;
}

/* ---------------------------------------------------------------------- */

int PythonImpl::execute_string(char *cmd)
{
  PyUtils::GIL lock;
  int err = PyRun_SimpleString(cmd);
  if (err) PyUtils::Print_Errors();
  return err;
}

/* ---------------------------------------------------------------------- */

int PythonImpl::execute_file(char *fname)
{
  FILE *fp = fopen(fname, "r");
  if (fp == nullptr) return -1;

  PyUtils::GIL lock;
  int err = PyRun_SimpleFile(fp, fname);
  if (err) PyUtils::Print_Errors();

  if (fp) fclose(fp);
  return err;
}

/* ------------------------------------------------------------------ */

void PythonImpl::deallocate(int i)
{
  delete[] pfuncs[i].itype;
  delete[] pfuncs[i].ivarflag;
  delete[] pfuncs[i].ivalue;
  delete[] pfuncs[i].dvalue;
  for (int j = 0; j < pfuncs[i].ninput; j++) delete[] pfuncs[i].svalue[j];
  delete[] pfuncs[i].svalue;
  delete[] pfuncs[i].internal_var;
  delete[] pfuncs[i].ovarname;
  delete[] pfuncs[i].longstr;
}

/* ------------------------------------------------------------------ */

bool PythonImpl::has_minimum_version(int major, int minor)
{
  return (PY_MAJOR_VERSION == major && PY_MINOR_VERSION >= minor) || (PY_MAJOR_VERSION > major);
}

/* ------------------------------------------------------------------ */

void PythonImpl::finalize()
{
  if (Py_IsInitialized()) Py_Finalize();
}
