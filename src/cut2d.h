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

#ifndef SPARTA_CUT2D_H
#define SPARTA_CUT2D_H

#include "pointers.h"
#include "my_vec.h"

namespace SPARTA_NS {

class Cut2d : protected Pointers {
 public:
  struct Cline {
    double x[2],y[2];   // coords of end points of line clipped to cell
    int line;           // index in list of lines that intersect this cell
  };

  struct Point {
    double x[2];        // coords of point
    int type;           // type of pt = ENTRY,EXIT,TWO,CORNER
                        // ENTRY/EXIT = only part of one Cline,
                        //   could also be a geometric corner pt
                        // TWO = part of two Clines
                        // CORNER = part of no Cline, is a geometric corner pt
    int next;           // index of next point when walking a loop
                        // set for ENTRY/TWO pts between ENTRY and EXIT
                        // set for EXIT/CORNER points around cell perimeter,
                        //   though may not be walked
    int line;           // original line (as stored by Cline) the pt starts,
                        //   only set for ENTRY and TWO pts
    int corner;         // 1,2,3,4 if x is a corner point, else 0
                        // could be ENTRY,EXIT,CORNER pt, but not a TWO pt
    int cprev,cnext;    // indices of pts in linked list around cell perimeter
    int side;           // which side of cell (0,1,2,3) pt is on
                        // only for ENTRY/EXIT/CORNER pts to make linked list
    double value;       // coord along the side
                        // only for ENTRY/EXIT/CORNER pts to make linked list
  };

  struct Loop {
    double area;        // area of loop
    int active;         // 1/0 if active or not
    int flag;           // INTERIOR (if all TWO points) or BORDER
    int n;              // # of points in loop
    int first;          // index of first point in loop
    int next;           // index of next loop in same PG, -1 if last loop
  };

  struct PG {
    double area;        // summed area (over loops) of PG
    int n;              // # of loops in PG
    int first;          // index of first loop in PG
  };

  MyVec<Cline> clines;  // list of Clines
  MyVec<Point> points;  // list of Points = Weiler/Atherton data structure
  MyVec<Loop> loops;    // list of loops in Points
  MyVec<PG> pgs;        // list of polygons = one or more loops

  Cut2d(class SPARTA *);
  ~Cut2d() {}
  int surf2grid(cellint, double *, double *, int *, int);
  int split(cellint, double *, double *, int, int *,
            double *&, int *, int *, int &, double *);
  void split_face(int, int, double *, double *);

 private:
  cellint id;            // ID of cell being worked on
  double *lo,*hi;        // opposite corner pts of cell
  int nsurf;             // # of surf elements in cell
  int *surfs;            // indices of surf elements in cell

  MyVec<double> areas;   // areas of each flow polygon found
  MyVec<int> used;       // 0/1 flag for each point when walking loops

  int build_clines();
  void weiler_build();
  void weiler_loops();
  void loop2pg();
  void create_surfmap(int *);
  int split_point(int *, double *);

  int cliptest(double *, double *);
  void clip(double *, double *, double *, double *);

  int ptflag(double *);
  int sameedge(double *, double *);
  int whichside(double *);

  void print_clines();
  void print_points();
  void print_loops();
};

}

#endif

/* ERROR/WARNING messages:

E: Bad grid of processors for create_grid

UNDOCUMENTED

E: Per-processor grid count is too big

UNDOCUMENTED

*/
