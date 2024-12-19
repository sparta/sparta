/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(ablate,FixAblate)

#else

#ifndef SPARTA_FIX_ABLATE_H
#define SPARTA_FIX_ABLATE_H

#include "fix.h"

namespace SPARTA_NS {

class FixAblate : public Fix {
 public:
  int igroup,dim;
  int nx,ny,nz;
  double cornerlo[3];
  double xyzsize[3];

  FixAblate(class SPARTA *, int, char **);
  ~FixAblate();
  int setmask();
  void init();
  void setup() {}
  void end_of_step();

  int pack_grid_one(int, char *, int);
  int unpack_grid_one(int, char *);
  void copy_grid_one(int, int);
  void reset_grid_count(int);
  void add_grid_one();
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

  void store_corners(int, int, int, double *, double *,
                     double **, double ***, int *, double, char *, int);
  double mindist;             // min fractional distance between any grid corner pt
                              //   and a generated tri vertex or line segment endpt
 protected:
  int me;
  int groupbit,which,argindex,icompute,ifix,ivariable,maxrandom;
  double scale;
  char *idsource;
  int storeflag;
  int multi_val_flag;
  int multi_dec_flag;
  int minmaxflag;
  int ncorner;
  int nmultiv;
  int sgroupbit;
  double thresh;
  double sum_delta;
  int ndelete;

  int nglocal;            // # of owned grid cells

  double **cvalues;       // corner point values
  double ***mvalues;      // corner multi values
  int *tvalues;           // per-cell type value
  int tvalues_flag;       // 1 if tvalues is defined (by ReadIsurf)

  double corner_inside_min;   // min allowed value for an inside corner point
  double corner_outside_max;  // max allowed value for an outside corner point

  int **ixyz;             // ix,iy,iz indices (1 to Nxyz) of my cells
                          // in 2d/3d ablate grid (iz = 1 for 2d)
  // DEBUG
  int **mcflags;

  double *celldelta;       // per-cell delta from compute or fix source
  double **cdelta;         // per-corner point deltas
  double **cdelta_ghost;   // ditto for my ghost cells communicated to me
  double ***mdelta;        // cdelta for multivalues
  double ***mdelta_ghost;  // ditto for my ghost cells (multivalues)
  double **nvert;          // number of vertices around each corner
  double **nvert_ghost;    // ditto for my ghost cells communicated to me

  int maxgrid;             // max size of per-cell vectors/arrays
  int maxghost;            // max size of cdelta_ghost

  int *proclist;
  cellint *locallist;
  int *numsend;
  int maxsend;

  double *sbuf;
  int maxbuf;

  double *vbuf;
  int maxvar;

  // new ablate with distributed decrement

  int refcorners[8];

  class MarchingSquares *ms;
  class MarchingCubes *mc;
  class RanKnuth *random;

  void process_args(int, char **);

  void create_surfs(int);

  void set_delta();
  void set_delta_random();
  void set_delta_uniform();

  void decrement();
  void decrement_multiv();
  void decrement_multid_inside();
  void decrement_multid_outside();
  void decrement_multiv_multid_outside();
  void decrement_multiv_multid_inside();

  void sync();
  void sync_multiv();
  void sync_multid_outside();
  void sync_multid_inside();
  void sync_multiv_multid_outside();
  void sync_multiv_multid_inside();

  void mark_corners_2d(int);
  void mark_corners_3d(int);
  int find_ninter();

  void epsilon_adjust();
  void epsilon_adjust_multiv();
  void push_lohi();
  void comm_neigh_corners(int);
  int walk_to_neigh(int, int, int, int);
  void grow_percell(int);
  void grow_send();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/
