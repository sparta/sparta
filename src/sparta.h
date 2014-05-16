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

#ifndef SPARTA_SPARTA_H
#define SPARTA_SPARTA_H

#include "stdio.h"

namespace SPARTA_NS {

class SPARTA {
 public:

  // fundamental SPARTA classes

  class Memory *memory;          // memory allocation functions
  class Error *error;            // error handling
  class Universe *universe;      // universe of processors
  class Input *input;            // input script processing
                                 // ptrs to top-level SPARTA-specific classes
  class Particle *particle;      // particles
  class Update *update;          // timestepper
  class Comm *comm;              // inter-processor communication
  class Domain *domain;          // simulation box
  class Modify *modify;          // fixes and computes
  class Grid *grid;              // volumetric grid cells
  class Surf *surf;              // surface elements
  class Collide *collide;        // collisions and chemistry
  class React *react;            // chemistry reactions
  class Output *output;          // stats/dump/restart
  class Timer *timer;            // CPU timing info

  MPI_Comm world;                // MPI communicator
  FILE *infile;                  // infile
  FILE *screen;                  // screen output
  FILE *logfile;                 // logfile

  // other top-level SPARTA classes and variables

  SPARTA(int, char **, MPI_Comm);
  ~SPARTA();
  void create();
  void init();
  void destroy();

  void print_styles();
};

}

#endif

/* ERROR/WARNING messages:

*/
