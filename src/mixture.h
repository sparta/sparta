/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_MIXTURE_H
#define SPARTA_MIXTURE_H

#include "pointers.h"

namespace SPARTA_NS {

class Mixture : protected Pointers {
 public:
  char *id;                   // ID of mixture
  int nspecies;               // # of species in mixture
  int *species;               // species indices in Particle species list

                              // global attributes
  double nrho;                // number density
  int nrho_flag;              // 1 if user set nrho
  double nrho_user;           // user value
  double vstream[3];          // stream velocity
  int vstream_flag;           // 1 if user set vstream
  double vstream_user[3];     // user value
  double temp_thermal;        // thermal temperature
  int temp_thermal_flag;      // 1 if user set thermal temp
  double temp_thermal_user;   // user value

                              // per-species attributes
  double *fraction;           // relative fraction of each species
  int *fraction_flag;         // 1 if user set fraction for a species
  double *fraction_user;      // user fractional value
  double *cummulative;        // cummulative fraction for each species

  int ngroup;                 // # of defined groups
  char **groups;              // group IDs
  int *mix2group;             // m2g[i] = group that mixture species I is in
  int *species2group;         // s2g[i] = group that Particle species I is in
                              // -1 if species I not in mixture

  double *vscale;             // pre-computed velocity scale factor

  Mixture(class SPARTA *, char *);
  ~Mixture();
  void command(int, char **);
  void init();
  void add_species_default(char *);
  
 private:
  int maxspecies,maxgroup;

  int *active;                // 1 if species is active for mixture command
  int nactive;                // # of active species for command
  int all_default;            // 1 if this is default mixture "all"
  int species_default;        // 1 if this is default mixture "species"

  void add_species(int, char **);
  void params(int, char **);
  void allocate();
  void delete_groups();
  void shrink_groups();
  void add_group(const char *);
  int find_group(const char *);
};

}

#endif
