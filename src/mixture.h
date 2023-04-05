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

#ifndef SPARTA_MIXTURE_H
#define SPARTA_MIXTURE_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class Mixture : protected Pointers {
 public:
  char *id;                   // ID of mixture
  int nspecies;               // # of species in mixture
  int *species;               // species[i] = particle species index of
                              //              mixture species I

  int ngroup;                 // # of defined groups
  char **groups;              // group IDs
  int *mix2group;             // m2g[i] = group that mixture species I is in

                              // global attributes
  double nrho;                // number density
  int nrho_flag;              // 1 if user set nrho
  double nrho_user;           // user value
  double vstream[3];          // stream velocity
  int vstream_flag;           // 1 if user set vstream
  double vstream_user[3];     // user value
  double temp_thermal;        // thermal temperature
  double temp_rot;            // rotational temperature
  double temp_vib;            // vibrational temperature
  int temp_thermal_flag;      // 1 if user set thermal temp
  int temp_rot_flag;          // 1 if user set rotational temp
  int temp_vib_flag;          // 1 if user set vibrational temp
  double temp_thermal_user;   // user value
  double temp_rot_user;       // user value
  double temp_vib_user;       // user value

                              // per-species attributes
  double *fraction;           // relative fraction of each species
  int *fraction_flag;         // 1 if user set fraction for a species
  double *fraction_user;      // user fractional value

  // set by init()

  double *cummulative;        // cummulative fraction for each species
  int *groupsize;             // # of species in each group
  int **groupspecies;         // list of particle species indices in each group
  int *species2group;         // s2g[i] = group that particle species I is in
                              // -1 if species I not in mixture
  int *species2species;       // s2s[i] = mixture species that
                              //   particle species I is
                              // -1 if species I not in mixture
  double *vscale;             // pre-computed velocity scale factor

  Mixture(class SPARTA *, char *);
  ~Mixture();
  void copy(Mixture *);
  void command(int, char **);
  void init();
  int init_fraction(int *, double *, double *, double *);
  void add_species_default(char *);
  int find_group(const char *);
  void write_restart(FILE *fp);
  void read_restart(FILE *fp);

 private:
  int maxspecies,maxgroup;
  int copyflag,copyarg;

  int activeflag;             // 1 if species are listed in mixture command
  int *active;                // flags for species listed in mixture command
  int all_default;            // 1 if this is default mixture "all"
  int species_default;        // 1 if this is default mixture "species"

  void add_species(int, char **);
  void params(int, char **);
  void allocate();
  void delete_groups();
  void shrink_groups();
  void add_group(const char *);
};

}

#endif

/* ERROR/WARNING messages:

E: Mixture ID must be alphanumeric or underscore characters

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Mixture %s fractions exceed 1.0

The sum of fractions must not be > 1.0.

E: Mixture species is not defined

One or more of the species ID is unknown.

E: Cannot add new species to mixture all or species

This is done automatically for these 2 mixtures when
each species is defined by the species command.

E: Mixture group ID must be alphanumeric or underscore characters

Self-explanatory.

E: Cannot use group keyword with mixture all or species

This is because the groups for these 2 mixtures are
pre-defined.

*/
