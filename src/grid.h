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
                              // set for all cells
    int id;                   // cell ID
    double lo[3],hi[3];       // opposite corner pts of cell
    int neigh[6];             // global indices of 6 neighbor cells
                              // XLO,XHI,YLO,YHI,ZLO,ZHI
                              // -1 if global boundary, including ZLO/ZHI in 2d
    int proc;                 // proc that owns this cell
    int nsurf;                // # of lines or triangles in this cell
    int nsplit;               // 1, unsplit leaf cell
                              // N > 1, split cell parent with N children
                              // N <= 0, split cell child, -N = parent

                              // set only for split parent cells
    double xsplit[3];         // coords of known point in parent split cell
    int xchild;               // index of child split cell that xsplit is in

                              // set only for owned cells
    int type;                 // CELLOUTSIDE,CELLINSIDE,CELLOVERLAP
    int plocal;               // local index of cell in myparent
    int clocal;               // local index of cell in mychild
    int count;                // # of particles in this cell
    int first;                // index of 1st particle in this cell
    double volume;            // flow volume of cell
                              // set for unsplit leaf and split cell child
  };

  OneCell *cells;             // global list of grid cells
  int ncell;                  // total # of original grid cells
  int nsplit;                 // total # of added split cells

  int *myparent;              // indices of parent grid cells I own
  int nparent;                // # of original grid cells I own
                              //   unsplit and split parents

  int *mychild;               // indices of child grid cells I own w/ particles
  int nchild;                 // # of child grid cells I own
                              //   unsplit parents and split children

  int **csurfs;      // indices of lines/tris in each cell
                     // ncell by cells->nsurf in size (ragged array)

  int **csplits;     // indices of which split cell each surf belongs to
                     // ncell by cells->nsurf in size (ragged array)
                     // NOTE: for now, a very inefficient data struct, 
                     //       only used for split parent cells

  int **cflags;      // SURFEXT,SURFINT,SURFOVERLAP for each cell corner point
                     // nlocal by 4/8 for 2d/3d
                     // corner pts ordered by x first, y next, z last

  Grid(class SPARTA *);
  ~Grid();
  void init();
  void add_cell(int, double *, double *, int *);
  void add_split_cell(int);
  int which_cell(double, double, double);
  void assign_stride(int);
  void assign_clump(int);
  void assign_block(int, int, int);
  void assign_random();
  void grow(int);
  bigint memory_usage();

 private:
  int maxcell;

  void assign_parent();
  void assign_child();
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
