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
  int sorted;               // 1 if particles are sorted by grid cell

  struct Species {          // info on each particle species
    char id[16];
    double molwt;
    double mass;
    double rotrel;
    double vibrel;
    double vibtemp;
    double specwt;
    double charge;
    int rotdof,vibdof;
    int internaldof;
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
    double evib;            // vibrational energy
    int flag;               // used for migration status
    double dtremain;        // portion of move timestep remaining
    double weight;          // particle or cell weight, if weighting enabled
  };

  struct OnePartRestart {
    int id;                 // particle ID
    int ispecies;           // particle species index
    cellint icell;          // cell ID the particle is in
    int nsplit;             // 1 for unsplit cell
                            // else neg of sub cell index (0 to Nsplit-1)
    double x[3];            // particle position
    double v[3];            // particle velocity
    double erot;            // rotational energy
    double evib;            // vibrational energy
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

  // extra custom vectors/arrays for per-particle data
  // ncustom > 0 if there are any extra arrays
  // custom attributes are created by various commands
  // these varaiables are public, others below are private

  int ncustom;              // # of custom attributes, some may be deleted
  int *etype;               // type = INT/DOUBLE of each attribute
  int *esize;               // size = 0 for vector, N for array columns
  int *ewhich;              // index into eivec,eiarray,edvec,edarray for data

  int **eivec;              // pointer to each integer vector
  int ***eiarray;           // pointer to each integer array
  double **edvec;           // pointer to each double vector
  double ***edarray;        // pointer to each double array

  // restart buffers, filled by read_restart

  int nlocal_restart;
  char *particle_restart;

  // methods

  Particle(class SPARTA *);
  ~Particle();
  void init();
  void compress_migrate(int, int *);
  void compress_rebalance();
  void compress_reactions(int, int *);
  void sort();
  void remove_all_from_cell(int);
  void grow(int);
  void grow_next();
  void pre_weight();
  void post_weight();

  int add_particle(int, int, int, double *, double *, double, double);
  int clone_particle(int);
  void add_species(int, char **);
  void add_mixture(int, char **);
  int find_species(char *);
  int find_mixture(char *);
  double erot(int, double, class RanPark *);
  double evib(int, double, class RanPark *);

  void write_restart_species(FILE *fp);
  void read_restart_species(FILE *fp);
  void write_restart_mixture(FILE *fp);
  void read_restart_mixture(FILE *fp);

  int size_restart();
  int pack_restart(char *);
  int unpack_restart(char *);

  int find_custom(char *);
  int add_custom(char *, int, int);
  void grow_custom(int, int, int);
  void remove_custom(int);
  void copy_custom(int, int);
  int sizeof_custom();
  void write_restart_custom(FILE *fp);
  void read_restart_custom(FILE *fp);
  void pack_custom(int, char *);
  void unpack_custom(char *, int);

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

  // extra custom vectors/arrays for per-particle data
  // ncustom > 0 if there are any extra arrays
  // these varaiables are private, others above are public

  char **ename;             // name of each attribute

  int ncustom_ivec;         // # of integer vector attributes
  int ncustom_iarray;       // # of integer array attributes
  int *icustom_ivec;        // index into ncustom for each integer vector
  int *icustom_iarray;      // index into ncustom for each integer array
  int *eicol;               // # of columns in each integer array (esize)

  int ncustom_dvec;         // # of double vector attributes
  int ncustom_darray;       // # of double array attributes
  int *icustom_dvec;        // index into ncustom for each double vector
  int *icustom_darray;      // index into ncustom for each double array
  int *edcol;               // # of columns in each double array (esize)

  int *custom_restart_flag; // flag on each custom vec/array read from restart
                            // used to delete them if not redefined in 
                            // restart script

  // private methods

  void read_species_file();
  int wordcount(char *, char **);
};

}

#endif

/* ERROR/WARNING messages:

E: Per-processor particle count is too big

No processor can have more particle than fit in a 32-bit integer,
approximately 2 billion.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Cannot open species file %s

Self-explanatory.

E: Invalid character in species ID

The only allowed characters are alphanumeric, an underscore, a plus
sign, or a minus sign.

E: Species ID is already defined

Species IDs must be unique.

E: Species ID does not appear in species file

Could not find the requested species in the specified file.

E: Incorrect line format in species file

Line read did not have expected number of fields.

E: Invalid species ID in species file

Species IDs are limited to 15 characters.

*/
