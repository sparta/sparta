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

#ifndef DSMC_DSMC_H
#define DSMC_DSMC_H

#include "stdio.h"

namespace DSMC_NS {

class DSMC {
 public:
                                 // ptrs to fundamental DSMC classes
  class Memory *memory;          // memory allocation functions
  class Error *error;            // error handling
  class Universe *universe;      // universe of processors
  class Input *input;            // input script processing
                                 // ptrs to top-level DSMC-specific classes
  class Particle *particle;      // particles
  class Update *update;          // timestepper
  class Comm *comm;              // inter-processor communication
  class Domain *domain;          // simulation box
  class Grid *grid;              // volumetric grid cells
  class Surf *surf;              // surface elements
  class Output *output;          // stats/dump/restart
  class Timer *timer;            // CPU timing info

  MPI_Comm world;                // MPI communicator
  FILE *infile;                  // infile
  FILE *screen;                  // screen output
  FILE *logfile;                 // logfile

  char *suffix;                  // suffix to add to input script style names
  int suffix_enable;             // 1 if suffix enabled, 0 if disabled
  class Cuda *cuda;              // CUDA accelerator class

  DSMC(int, char **, MPI_Comm);
  ~DSMC();
  void create();
  void init();
  void destroy();

  void print_styles();
};

}

#endif
