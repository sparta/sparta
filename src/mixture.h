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
  double nrho;
  double nrho_user;
  double *fraction;
  double **vstream;
  double *temp_thermal;
  double *fraction_user;
  double **vstream_user;
  double *temp_thermal_user;
  int nrho_flag;
  int *fraction_flag;
  int *vstream_flag;
  int *temp_thermal_flag;
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
