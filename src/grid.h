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

#ifndef SPARTA_GRID_H
#define SPARTA_GRID_H

#include "pointers.h"

namespace SPARTA_NS {

class Grid : protected Pointers {
 public:
  int grid_exist;

  int nx,ny,nz;
  double xdelta,ydelta,zdelta;
  double xdeltainv,ydeltainv,zdeltainv;

  struct OneCell {
    int id;
    int type;                 // SURFEXTERIOR,SURFINTERIOR,SURFOVERLAP (owned)
    double lo[3],hi[3];       // opposite corner pts of cell
    int neigh[6];             // global indices of 6 neighbor cells
                              // XLO,XHI,YLO,YHI,ZLO,ZHI
                              // -1 if global boundary, including ZLO/ZHI in 2d
    int proc;                 // proc that owns this cell
    int local;                // local index of cell (owned)
    int count;                // # of particles in this cell (owned)
    int first;                // index of 1st particle in this cell (owned)
    int nsurf;                // # of lines or triangles in this cell
    double volume;            // flow volume of cell (owned)
                              // set for unsplit leaf or split cell child
    int nsplit;               // 1, unsplit leaf cell
                              // N > 1, split cell parent with N children
                              // N <= 0, split cell child, -N = index of parent
  };

  OneCell *cells;             // global list of grid cells
  int ncell;                  // total # of grid cells

  int *mycells;               // indices of grid cells I own
  int nlocal;                 // # of grid cells I own
  
  int **csurfs;      // indices of lines/tris in each cell
                     // ncell by cells->nsurf in size (ragged array)

  int **csplits;     // for split cells only
                     // indices of which split cell each surf belongs to
                     // ncell by cells->nsurf in size (ragged array)

  int **cflags;      // SURFEXT,SURFINT,SURFOVERLAP for each cell corner point
                     // nlocal by 4/8 for 2d/3d
                     // corner pts ordered by x first, y next, z last

  Grid(class SPARTA *);
  ~Grid();
  void init();
  void add_cell(int, double *, double *, int *);
  int which_cell(double, double, double);
  void assign_stride(int);
  void assign_block(int, int, int);
  void assign_random();
  void grow(int);
  bigint memory_usage();

 private:
  int maxcell;

  void assign_mine();
  void procs2grid(int &, int &, int &);
  void surf2grid();
  void grid_inout();
  void grid_inout2();
  int flood(int, int, int);
  void grid_check();
  void surf2grid_stats();
  void flow_stats();
  double flow_volume();
};

}

#endif

/* ERROR/WARNING messages:

E: Bad grid of processors for create_grid

UNDOCUMENTED

E: Per-processor grid count is too big

UNDOCUMENTED

*/
