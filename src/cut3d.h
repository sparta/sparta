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

#ifndef SPARTA_CUT3D_H
#define SPARTA_CUT3D_H

#include "pointers.h"
#include "cut2d.h"
#include "my_vec.h"

namespace SPARTA_NS {

class Cut3d : protected Pointers {
 public:
  int ntiny,nshrink;

  Cut3d(class SPARTA *);
  ~Cut3d();
  int surf2grid(cellint, double *, double *, surfint *, int);
  int surf2grid_list(cellint, double *, double *, int, surfint *,
                     surfint *, int);
  int surf2grid_one(double *, double *, double *, double *, double *);

  int split(cellint, double *, double *, int, surfint *,
            double *&, int *, int *, int &, double *);

  int clip_external(double *, double *, double *,
                    double *, double *, double *);
  int sameface(double *, double *, double *);
  int sameface_external(double *, double *, double *, double *, double *);

 private:
  int implicit;

  cellint id;            // ID of cell being worked on
  double lo[3],hi[3];    // opposite corner pts of cell being worked on

  int nsurf;             // # of surf elements in cell
  surfint *surfs;        // indices of surf elements in cell, caller owns

  int grazecount;        // count of tris that graze cell surf w/ outward norm
  int touchcount;        // count of tris that only touch cell surf
  int touchmark;         // corner marking inferred by touching tris
  double epsilon;        // epsilon size for this cell

  double **path1,**path2;

  MyVec<double> vols;    // vols of each flow polyhedron found

  int empty;

  struct Vertex {
    int active;      // 1/0 if active or not
    int style;       // CTRI or CTRIFACE or FACEPGON or FACE
    int label;       // index in list of tris that intersect this cell
                     //   for CTRI or CTRIFACE
                     // face index (0-5) for FACEPGON or FACE
    int next;        // index of next vertex when walking a loop
    int nedge;       // # of edges in this vertex
    int first;       // first edge in vertex
    int dirfirst;    // dir of first edge in vertex
    int last;        // last edge in vertex
    int dirlast;     // dir of last edge in vertex
    double volume;   // volume of vertex projected against lower z face of cell
    double *norm;    // ptr to norm of tri, NULL for other styles
  };

  struct Edge {
    double p1[3],p2[3];  // 2 points in edge
    int active;          // 1/0 if active or not
    int style;           // CTRI or CTRIFACE or FACEPGON or FACE
    int clipped;         // 1/0 if already clipped during face iteration
    int nvert;           // flag for verts containing this edge
                         // 0 = no verts
                         // 1 = just 1 vert in forward dir
                         // 2 = just 1 vert in reverse dir
                         // 3 = 2 verts in both dirs
                         // all vecs are [0] in forward dir, [1] in reverse dir
    int verts[2];        // index of vertices containing this edge, -1 if not
    int next[2];         // index of next edge for each vertex, -1 for end
    int dirnext[2];      // next edge for each vertex is forward/reverse (0,1)
    int prev[2];         // index of prev edge for each vertex, -1 for start
    int dirprev[2];      // prev edge for each vertex is forward/reverse (0,1)
  };

  struct Loop {
    double volume;        // volume of loop
    int flag;             // INTERIOR (if all CTRI vertices) or BORDER
    int n;                // # of vertices in loop
    int first;            // index of first vertex in loop
    int next;             // index of next loop in same PH, -1 if last loop
  };

  struct PH {
    double volume;
    int n;
    int first;
  };

  MyVec<Vertex> verts;       // list of vertices in BPG
  MyVec<Edge> edges;         // list of edges in BPG
  MyVec<Loop> loops;         // list of loops of vertices = polyhedra
  MyVec<PH> phs;             // list of polyhedrons = one or more loops

  MyVec<int> facelist[6];    // list of edges on each cell face
  MyVec<int> used;           // 0/1 flag for each vertex when walking loops
  MyVec<int> stack;          // list of vertices to check when walking loops

  class Cut2d *cut2d;

  // methods

  int clip(double *, double *, double *);
  int split_try(cellint, int, surfint *,
                double *&, int *, int *, int &, double *, int &);
  void split_error(int);

  int add_tris();
  void clip_tris();
  void clip_adjust();
  void ctri_volume();
  int edge2face();
  void edge2clines(int);
  int add_face_pgons(int);
  int add_face(int, double *, double *);
  void remove_faces();
  int check();
  void walk();
  int loop2ph();
  void create_surfmap(int *);
  int split_point_explicit(int *, double *, int &);
  int split_point_implicit(int *, double *, int &);

  void edge_insert(int, int, int, int, int, int, int);
  void edge_remove(Edge *);
  void edge_remove(Edge *, int);
  void vertex_remove(Vertex *);
  int grazing(Vertex *);
  int which_faces(double *, double *, int *);
  void face_from_cell(int, double *, double *);
  void compress2d(int, double *, double *);
  void expand2d(int, double, double *, double *);

  int findedge(double *, double *, int, int &);
  void between(double *, double *, int, double, double *);
  int samepoint(double *, double *);
  int corner(double *);
  int on_faces(double *, int *);
  void move_to_faces(double *);
  int ptflag(double *);

  void failed_cell();
  void print_bpg(const char *);
  void print_loops();
};

}

#endif

/* ERROR/WARNING messages:

E: Singlet BPG edge not on cell face

This is an error when calculating how a 3d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: BPG edge on more than 2 faces

This is an error when calculating how a 3d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: Vertex has less than 3 edges

This is an error when calculating how a 3d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: Vertex contains invalid edge

This is an error when calculating how a 3d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: Vertex contains edge that doesn't point to it

This is an error when calculating how a 3d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: Vertex contains duplicate edge

This is an error when calculating how a 3d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: Vertex pointers to last edge are invalid

This is an error when calculating how a 3d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: Edge not part of 2 vertices

This is an error when calculating how a 3d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: Edge part of same vertex twice

This is an error when calculating how a 3d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: Edge part of invalid vertex

This is an error when calculating how a 3d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: No positive volumes in cell

This is an error when calculating how a 3d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: More than one positive volume with a negative volume

SPARTA cannot determine which positive volume the negative volume is
inside of, if a cell is so large that it includes both positive and
negative volumes.

E: Single volume is negative, inverse donut

An inverse donut is a surface with a flow region interior to the donut
hole and also exterior to the entire donut.  This means the flow
regions are disconnected.  SPARTA cannot correctly compute the flow
volume of this kind of object.

E: Could not find split point in split cell

This is an error when calculating how a grid cell is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

E: Found edge in same direction

This is an error when calculating how a 3d grid is cut or split by
surface elements.  It should not normally occur.  Please report the
issue to the SPARTA developers.

*/
