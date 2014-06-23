/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov
   Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_SURF_H
#define SPARTA_SURF_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class Surf : protected Pointers {
 public:
  int exist;                // 1 if any surfaces are defined, else 0
  double bblo[3],bbhi[3];   // bounding box around surfs

  int nid;                  // # of unique surface IDs
  char **ids;               // list of IDs
  int *idlo,*idhi;          // range of consecutive surfs assigned to each ID

  struct Point {
    double x[3];
  };

  struct Line {
    int isc;                // index of surface collision model it belongs to
                            // -1 if unassigned
    int p1,p2;              // indices of points in line segment
                            // rhand rule: Z x (p2-p1) = outward normal
    double norm[3];         // outward normal to line segment
  };

  struct Tri {
    int isc;                // index of surface collision model it belongs to
                            // -1 if unassigned
    int p1,p2,p3;           // indices of points in triangle
                            // rhand rule: (p2-p1) x (p3-p1) = outward normal
    double norm[3];         // outward normal to triangle
  };

  Point *pts;               // global list of points
  Line *lines;              // global list of lines
  Tri *tris;                // global list of tris
  int npoint,nline,ntri;    // number of each

  int *mysurfs;             // indices of surf elements I own
  int nlocal;               // # of surf elements I own

  class SurfCollide **sc;      // list of surface collision models
  int nsc;                     // # of surface collision models
  int maxsc;                   // max # of models in sc

  Surf(class SPARTA *);
  ~Surf();
  void modify_params(int, char **);
  void init();
  int nelement();
  void setup_surf();

  void compute_line_normal(int, int);
  void compute_tri_normal(int, int);
  void quad_corner_point(int, double *, double *, double *);
  void hex_corner_point(int, double *, double *, double *);
  double line_size(int);
  double tri_size(int, double &);

  int add_surf(const char *);
  int find_surf(const char *);
  void add_collide(int, char **);
  int find_collide(const char *);

  void collate_vec(int, int *, double *, int, double *, int, int);
  void collate_array(int, int, int *, double **, double **);

  void write_restart(FILE *);
  void read_restart(FILE *);
  bigint memory_usage();
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Could not find surf_modify surf-ID

Self-explanatory.

E: Could not find surf_modify sc-ID

Self-explanatory.

E: %d surface elements not assigned to a collision model

All surface elements must be assigned to a surface collision model via
the surf_modify command before a simulation is perforemd.

E: Reuse of surf_collide ID

A surface collision model ID cannot be used more than once.

E: Invalid surf_collide style

Self-explanatory.

*/
