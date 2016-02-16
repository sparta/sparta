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
  int surf_collision_check; // flag for whether init() check is required
                            // for assign of collision models to surfs

  double bblo[3],bbhi[3];   // bounding box around surfs
  int tally_comm;           // style of comm for surf tallies

  int nreact_one;           // surface reactions in current step
  bigint nreact_running;    // running count of surface reactions

  int ngroup;               // # of defined groups
  char **gnames;            // name of each group
  int *bitmask;             // one-bit mask for each group
  int *inversemask;         // inverse mask for each group

  struct Point {
    double x[3];
  };

  struct Line {
    int type,mask;          // type and mask of the element
    int isc,isr;            // index of surface collision and reaction models
                            // -1 if unassigned
    int p1,p2;              // indices of points in line segment
                            // rhand rule: Z x (p2-p1) = outward normal
    double norm[3];         // outward normal to line segment
  };

  struct Tri {
    int type,mask;          // type and mask of the element
    int isc,isr;            // index of surface collision and reaction models
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

  int nsc,nsr;              // # of surface collision and reaction models
  class SurfCollide **sc;   // list of surface collision models
  class SurfReact **sr;     // list of surface reaction models

  int pushflag;             // set to 1 to push surf pts near grid cell faces
  double pushlo,pushhi;     // lo/hi ranges to push on
  double pushvalue;         // new position to push to

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
  double axi_line_size(int);
  double tri_size(int, double &);

  void check_watertight_2d(int, int);
  void check_watertight_3d(int, int);
  void check_point_inside(int, int);

  void add_collide(int, char **);
  int find_collide(const char *);
  void add_react(int, char **);
  int find_react(const char *);

  void group(int, char **);
  int add_group(const char *);
  int find_group(const char *);

  void collate_vector(int, int *, double *, int, double *);
  void collate_array(int, int, int *, double **, double **);

  void write_restart(FILE *);
  void read_restart(FILE *);
  bigint memory_usage();

 private:
  int maxsc;                // max # of models in sc
  int maxsr;                // max # of models in sr

  void collate_vector_allreduce(int, int *, double *, int, double *);
  void collate_vector_irregular(int, int *, double *, int, double *);
  void collate_array_allreduce(int, int, int *, double **, double **);
  void collate_array_irregular(int, int, int *, double **, double **);
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
