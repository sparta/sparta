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

  int nx,ny,nz;
  double xdelta,ydelta,zdelta;
  double xdeltainv,ydeltainv,zdeltainv;

  struct OneCell {
    int id;
    double lo[3],hi[3];       // opposite corner pts of cell
    int neigh[6];             // global indices of 6 neighbor cells
                              // -1 if global boundary
    int proc;                 // proc that owns this cell
    int local;                // local index of cell if I own it
    int first,count;
    double volume;            // volume of cell
  };

  OneCell *cells;             // global list of grid cells
  int ncell;                  // total # of grid cells

  int *mycells;               // indices of grid cells I own
  int nlocal;                 // # of grid cells I own
  
  Grid(class DSMC *);
  ~Grid();
  void init() {}
  void add_cell(int, double *, double *, int *);
  void setup_grid();
  int which_cell(double, double, double);
  void assign_stride(int);
  void assign_block(int, int, int);
  void assign_random(int);
  void grow(int);
  bigint memory_usage();

 private:
  int maxcell;

  void procs2grid(int &, int &, int &);
};

}

#endif

/* ERROR/WARNING messages:

E: Bad grid of processors for create_grid

UNDOCUMENTED

E: Per-processor grid count is too big

UNDOCUMENTED

*/
