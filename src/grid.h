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

#ifndef DSMC_GRID_H
#define DSMC_GRID_H

#include "pointers.h"

namespace DSMC_NS {

class Grid : protected Pointers {
 public:
  int grid_exist;

  struct OneCell {
    int id;
    double lo[3],hi[3];
    int neigh[6];
    int proc;
    int nparticles;
    int first;
  };

  OneCell *cells;
  int ncell;
  int nlocal;
  
  Grid(class DSMC *);
  ~Grid();
  void init() {}
  void create(int, char **);
  int which_cell(double, double, double);
  bigint memory_usage();

 private:
  int nx,ny,nz;
  double xdelta,ydelta,zdelta;
  double xdeltainv,ydeltainv,zdeltainv;
  int bstyle;
  int user_procgrid[3];
  int order,seed;

  int procgrid[3];

  void assign_stride();
  void assign_block();
  void assign_random();
  void procs2grid();
};

}

#endif
