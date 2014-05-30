/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_PARTICLE_H
#define SPARTA_PARTICLE_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class Particle : protected Pointers {
 public:
  int exist;                // 1 if particles exist

  struct Species {          // info on each particle species
    char id[16];
    double molwt;
    double mass;
    int rotdof;
    int rotrel;
    int vibdof;
    int vibrel;
    double vibtemp;
    double specwt;
    double charge;
  };

  Species *species;         // list of particle species info
  int nspecies;             // # of defined species

  class Mixture **mixture;
  int nmixture;
  int maxmixture;

  struct OnePart {
    int id;                 // particle ID
    int ispecies;           // particle species index
    int icell;              // which local Grid::cells the particle is in
    double x[3];            // particle position
    double v[3];            // particle velocity
    double erot;            // rotational energy
    int ivib;               // vibrational mode
    int flag;               // used for migration status
    double dtremain;        // portion of move timestep remaining
    double weight;          // particle or cell weight, if weighting enabled
  };

  struct OnePartRestart {
    int id;                 // particle ID
    int ispecies;           // particle species index
    cellint icell;          // cell ID the particle is in
    double x[3];            // particle position
    double v[3];            // particle velocity
    double erot;            // rotational energy
    int ivib;               // vibrational mode
  };

  bigint nglobal;           // global # of particles
  int nlocal;               // # of particles I own
  int maxlocal;             // max # particles list can hold
  OnePart *particles;       // list of particles I own

  // currently stored in grid.h for every cell, whether I own it or not
  // not sure why storing it here is slower

  //int *cellcount;           // count of particles in each grid cell I own
  //int *first;               // index of first particle in each grid cell

  int *next;                // index of next particle in each grid cell

  Particle(class SPARTA *);
  ~Particle();
  void init();
  void compress(int, int *);
  void compress();
  void sort();
  void grow(int);
  void pre_weight();
  void post_weight();

  int add_particle(int, int, int, double *, double *, double, int);
  int clone_particle(int);
  void add_species(int, char **);
  void add_mixture(int, char **);
  int find_species(char *);
  int find_mixture(char *);
  double erot(int, class RanPark *);
  int evib(int, class RanPark *);

  void write_restart(FILE *fp);
  void read_restart(FILE *fp);
  int size_restart();
  int pack_restart(char *);
  int unpack_restart(char *);

  bigint memory_usage();

 private:
  int me;
  int maxgrid;              // max # of indices first can hold
  int maxsort;              // max # of particles next can hold
  int maxspecies;           // max size of species list

  Species *filespecies;     // list of species read from file
  int nfilespecies;         // # of species read from file
  int maxfilespecies;       // max size of filespecies list
  FILE *fp;                 // species file pointer

  class RanPark *wrandom;   // RNG for particle weighting

  void read_species_file();
  int wordcount(char *, char **);
};

}

#endif

/* ERROR/WARNING messages:

E: Per-processor grid count is too big

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Cannot open species file %s

UNDOCUMENTED

E: Species ID is already defined

UNDOCUMENTED

E: Species ID does not appear in species file

UNDOCUMENTED

E: Incorrect line format in species file

UNDOCUMENTED

E: Invalid species ID in species file

UNDOCUMENTED

*/
