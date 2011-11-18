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

#ifndef DSMC_PARTICLE_H
#define DSMC_PARTICLE_H

#include "stdio.h"
#include "pointers.h"

namespace DSMC_NS {

class Particle : protected Pointers {
 public:
  struct Species {          // info on each particle species
    char id[16];
    double molwt;
    double mass;
    double diam;
    int rotdof;
    int rotrel;
    int vibdof;
    int vibrel;
    double vibtemp;
    double specwt;
    double charge;
    double omega;
    double tref;
    double alpha;
  };

  Species *species;         // list of particle species info
  int nspecies;             // # of defined species

  struct OnePart {
    int id,type;            // particle ID, type
    int icell;              // local grid cell the particle is in (0 to N-1)
    double x[3];            // coords of particle
    double v[3];            // velocity of particle
  };

  bigint nglobal;           // global # of particles
  int nlocal;               // # of particles I own
  int maxlocal;             // max # particles list can hold
  OnePart *particles;       // list of particles I own

  int *cellcount;           // count of particles in each grid cell I own
  int *first;               // index of first particle in each grid cell
  int *next;                // index of next particle in each grid cell

  Particle(class DSMC *);
  ~Particle();
  void init();
  void add_particle(int, int, int, double, double, double);
  void add_species(int, char **);
  void compress(int, int *);
  void sort();
  void grow(int);
  bigint memory_usage();

 private:
  int maxgrid;              // max # of indices first can hold
  int maxsort;              // max # of particles next can hold
  int maxspecies;           // max size of species list

  Species *filespecies;     // list of species read from file
  int nfilespecies;         // # of species read from file
  int maxfilespecies;       // max size of filespecies list
  FILE *fp;                 // species file pointer

  void read_species_file();
  int wordcount(char *, char **);
};

}

#endif
