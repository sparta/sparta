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

#ifndef DSMC_MIXTURE_H
#define DSMC_MIXTURE_H

#include "pointers.h"

namespace DSMC_NS {

class Mixture : protected Pointers {
 public:
  char *id;
  int nspecies;
  int *species;
  double *fraction;
  double *fraction_user;
  double nrho;
  double nrho_user;
  double vstream[3];
  double vstream_user[3];
  double temp_thermal;
  double temp_thermal_user;
  int *fraction_flag;
  int nrho_flag;
  int vstream_flag;
  int temp_thermal_flag;

  int allspecies;             // 1 if mixture contains all species in model
  int ngroups;                // # of defined groups
  int *species2group;         // s2g[i] = map of species I (1 to Nsp) to
                              //   group (0 to Ngroups-1),
                              //   -1 if species not in mixture
  char **groups;              // group IDs

  double *cummulative;
  double *vscale;
  int *active;

  Mixture(class DSMC *, char *);
  ~Mixture();
  void init();
  void add_species(int, char **);
  void params(int, char **);

 private:
  int maxspecies;

  void allocate(int);
};

}

#endif
