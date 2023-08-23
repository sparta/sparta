/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
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
                        // TWO if part of two Clines
                        // ENTRY/EXIT if only part of one Cline
                        // CORNER = part of no Cline, geometric corner pt
    int next;           // index of next point when walking a flow area loop
                        // set for ENTRY/TWO pts between ENTRY and EXIT
                        // set for EXIT/CORNER points around cell perimeter,
                        //   though may not be walked
    int line;           // index of line this pt starts in intersecting line list
                        //   only set for ENTRY and TWO pts
    int corner;         // 1,2,3,4 if pt is geometrically a corner point, else 0
                        // could be ENTRY,EXIT,CORNER pt, but not a TWO pt
    int cprev,cnext;    // indices of pts in linked list around cell perimeter
                        //   walking in counter-clockwise direction
    int side;           // which side of cell (0,1,2,3) pt is on
                        // only for ENTRY/EXIT/CORNER pts to make linked list
    double value;       // coord along the side
                        // only for ENTRY/EXIT/CORNER pts to make linked list
  };

  struct Loop {
    double area;        // area of loop
    int active;         // 1/0 if active or not
    int flag;           // INTERIOR, if all TWO pts (TWO pt can be on border)
                        // BORDER, if all CORNER pts
                        // INTBORD, if otherwise
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

  Cut2d(class SPARTA *, int);
  ~Cut2d() {}
  int surf2grid(cellint, double *, double *, surfint *, int);
  int surf2grid_list(cellint, double *, double *, int, surfint *,
                     surfint *, int);
  int surf2grid_one(double *, double *, double *, double *);

  int split(cellint, double *, double *, int, surfint *,
            double *&, int *, int *, int &, double *);
  int split_face(int, int, double *, double *);

  int clip_external(double *, double *, double *, double *, double *);
  int sameedge(double *, double *);
  int sameedge_external(double *, double *, double *, double *);

 private:
  int axisymmetric;
  int implicit;

  cellint id;            // ID of cell being worked on
  double *lo,*hi;        // opposite corner pts of cell
  int nsurf;             // # of surf elements in cell
  surfint *surfs;        // indices of surf elements in cell

  int grazecount;        // count of lines that graze cell surf w/ outward norm
  int touchcount;        // count of line that only touch cell surf
  int touchmark;         // corner marking inferred by touching lines

  MyVec<double> areas;   // areas of each flow polygon found
  MyVec<int> used;       // 0/1 flag for each point when walking loops

  // methods

  void build_clines();
  int weiler_build();
  void weiler_loops();
  int loop2pg();
  void create_surfmap(int *);
  int split_point_explicit(int *, double *, int &);
  int split_point_implicit(int *, double *, int &);

  int cliptest(double *, double *);
  void clip(double *, double *, double *, double *);

  int ptflag(double *);
  int whichside(double *);

  void failed_cell();
  void print_clines();
  void print_points();
  void print_loops();
};

}

#endif

/* ERROR/WARNING messages:

E: Point appears first in more than one CLINE

This is an error when calculating how a 2d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.
SPARTA developers.

E: Point appears last in more than one CLINE

This is an error when calculating how a 2d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: Singlet CLINES point not on cell border

This is an error when calculating how a 2d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: No positive areas in cell

This is an error when calculating how a 2d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: More than one positive area with a negative area

SPARTA cannot determine which positive area the negative area is
inside of, if a cell is so large that it includes both positive and
negative areas.

E: Single area is negative, inverse donut

An inverse donut is a surface with a flow region interior to the donut
hole and also exterior to the entire donut.  This means the flow
regions are disconnected.  SPARTA cannot correctly compute the flow
area of this kind of object.

E: Could not find split point in split cell

This is an error when calculating how a grid cell is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

*/
