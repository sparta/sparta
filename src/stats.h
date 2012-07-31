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

#ifndef DSMC_STATS_H
#define DSMC_STATS_H

#include "pointers.h"

namespace DSMC_NS {

class Stats : protected Pointers {
 public:
  Stats(class DSMC *);
  ~Stats();
  void init();
  void modify_params(int, char **);
  void set_fields(int, char **);
  void header();
  void compute(int);
  int evaluate_keyword(char *, double *);

 private:
  char *line;
  char **keyword;
  int *vtype;

  int nfield;
  int me;

  char **format,**format_user;
  char *format_float_one_def;
  char *format_int_one_def;
  char *format_float_user,*format_int_user,*format_bigint_user;
  char format_bigint_one_def[8];

  int firststep;
  int flushflag,lineflag;

  double last_tpcpu,last_spcpu;
  double last_time;
  bigint last_step;

  bigint nmol;
                         // data used by routines that compute single values
  int ivalue;            // integer value to print
  double dvalue;         // double value to print
  bigint bivalue;        // big integer value to print
  int ifield;            // which field in thermo output is being computed
  int *field2index;      // which compute,fix,variable calcs this field
  int *argindex1;        // indices into compute,fix scalar,vector
  int *argindex2;

  int ncompute;                // # of Compute objects called by thermo
  char **id_compute;           // their IDs
  int *compute_which;          // 0/1/2 if should call scalar,vector,array
  class Compute **computes;    // list of ptrs to the Compute objects

  int nfix;                    // # of Fix objects called by thermo
  char **id_fix;               // their IDs
  class Fix **fixes;           // list of ptrs to the Fix objects

  int nvariable;               // # of variables evaulated by thermo
  char **id_variable;          // list of variable names
  int *variables;              // list of Variable indices

  // private methods

  void allocate();
  void deallocate();

  int add_compute(const char *, int);
  int add_fix(const char *);
  int add_variable(const char *);

  typedef void (Stats::*FnPtr)();
  void addfield(const char *, FnPtr, int);
  FnPtr *vfunc;                // list of ptrs to functions

  void compute_compute();      // functions that compute a single value
  void compute_fix();          // via calls to  Compute,Fix,Variable classes
  void compute_variable();

  // functions that compute a single value
  // customize a new keyword by adding a method prototype

  void compute_step();
  void compute_elapsed();
  void compute_dt();
  void compute_cpu();
  void compute_tpcpu();
  void compute_spcpu();

  void compute_nmol();
  void compute_ntouch();
  void compute_ncomm();
  void compute_nbound();
  void compute_nexit();
  void compute_nscoll();
  void compute_nscheck();
  void compute_ncoll();
  void compute_natt();

  void compute_nmolave();
  void compute_ntouchave();
  void compute_ncommave();
  void compute_nboundave();
  void compute_nexitave();
  void compute_nscollave();
  void compute_nscheckave();
  void compute_ncollave();
  void compute_nattave();

  void compute_vol();
  void compute_lx();
  void compute_ly();
  void compute_lz();

  void compute_xlo();
  void compute_xhi();
  void compute_ylo();
  void compute_yhi();
  void compute_zlo();
  void compute_zhi();
};

}

#endif
