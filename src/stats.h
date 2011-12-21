/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
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
  void header();
  void compute(int);
  void set_fields(int, char **);
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

  bigint npart;
                         // data used by routines that compute single values
  int ivalue;            // integer value to print
  double dvalue;         // double value to print
  bigint bivalue;        // big integer value to print
  int ifield;            // which field in thermo output is being computed
  int *field2index;      // which compute,fix,variable calcs this field
  int *argindex1;        // indices into compute,fix scalar,vector
  int *argindex2;

                         // data for keyword-specific Compute objects
                         // index = where they are in computes list
                         // id = ID of Compute objects
                         // Compute * = ptrs to the Compute objects

  // private methods

  void allocate();
  void deallocate();

  typedef void (Stats::*FnPtr)();
  void addfield(const char *, FnPtr, int);
  FnPtr *vfunc;                // list of ptrs to functions

  // functions that compute a single value
  // customize a new keyword by adding a method prototype

  void compute_step();
  void compute_elapsed();
  void compute_dt();
  void compute_cpu();
  void compute_tpcpu();
  void compute_spcpu();

  void compute_npart();
  void compute_temp();

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
