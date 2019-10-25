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
#include "hash3.h"
#include "hashlittle.h"

namespace SPARTA_NS {

class Surf : protected Pointers {
 public:
  int exist;                // 1 if any surfaces are defined, else 0
  int implicit;             // 1 = implicit surfs, 0 = explicit surfs
  int distributed;          // 1 = surfs spread across procs (exp or impl)
                            // 0 = each proc owns all
  //int dynamic;              // 1 = implicit surfs change, 0 = static
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

  // surf data structs
  // explicit, all: each proc owns all surfs
  //   nlocal = Nsurf, nghost = 0, tris = all surfs
  //   nown = Nsurf/P, mytris = NULL
  // explicit, distributed: each proc owns N/P surfs
  //   nlocal/nghost = surfs in owned/ghost grid cells, tris = nloc+ngh surfs
  //   nown = Nsurf/P, mytris = surfs I uniquely own
  // implicit, distributed: each proc owns surfs in its owned grid cells
  //   nlocal = surfs in owned grid cells, nghost = surfs in ghost grid cells
  //   nown = nlocal, mytris = NULL

  bigint nsurf;             // total # of surf elements, lines or tris

  struct Line {
    surfint id;             // unique ID for explicit surf
                            // cell ID for implicit surf
    int type,mask;          // type and mask of the surf
    int isc,isr;            // index of surface collision and reaction models
                            // -1 if unassigned
    double p1[3],p2[3];     // end points of line segment
                            // rhand rule: Z x (p2-p1) = outward normal
    double norm[3];         // outward normal to line segment
  };

  struct Tri {
    surfint id;             // unique ID for explicit surf
                            // cell ID for implicit surf
    int type,mask;          // type and mask of the surf
    int isc,isr;            // index of surface collision and reaction models
                            // -1 if unassigned
    double p1[3],p2[3],p3[3];  // corner points of triangle
                            // rhand rule: (p2-p1) x (p3-p1) = outward normal
    double norm[3];         // outward normal to triangle
  };

  Line *lines;              // list of lines for surface collisions
  Tri *tris;                // list of tris for surface collisions
  int nlocal;               // # of lines or tris
                            // explicit, all: nlocal = nsurf
                            // explicit, distributed: 
                            //   surfs which overlap my owned grid cells
                            // implicit: surfs within my owned grid cells
  int nghost;               // # of ghost surfs I store for collisions
                            // explicit, all: nghost = 0
                            // explicit, distributed: 
                            //   surfs which overlap my ghost grid cells
                            // implicit: surfs within my ghost grid cells
  int nmax;                 // max length of lines/tris vecs

  Line *mylines;            // list of lines assigned uniquely to me
                            //   only for explicit, distributed
  Tri *mytris;              // list of tris assigned uniquely to me
                            //   only for explicit, distributed
  int nown;                 // # of lines or tris I own uniquely
  int maxown;               // max length of owned lines/tris vecs

  Line *tmplines;           // list of temporary lines, filled by ReadSurf
  Tri *tmptris;             // list of temporary tris, filled by ReadSurf
  int ntmp;                 // # of temporary surfs
  int nmaxtmp;              // max size of tmplines/tmptris

  int nsc,nsr;              // # of surface collision and reaction models
  class SurfCollide **sc;   // list of surface collision models
  class SurfReact **sr;     // list of surface reaction models

  int pushflag;             // set to 1 to push surf pts near grid cell faces
  double pushlo,pushhi;     // lo/hi ranges to push on
  double pushvalue;         // new position to push to

#include "hash_options.h"

#ifdef SPARTA_MAP
  typedef std::map<surfint,int> MySurfHash;
  typedef std::map<OnePoint2d,int> MyHashPoint;
  typedef std::map<OnePoint2d,int>::iterator MyPointIt;
  typedef std::map<TwoPoint3d,int> MyHash2Point;
  typedef std::map<TwoPoint3d,int>::iterator My2PointIt;
#elif defined SPARTA_UNORDERED_MAP
  typedef std::unordered_map<surfint,int> MySurfHash;
  typedef std::unordered_map<OnePoint2d,int,OnePoint2dHash> MyHashPoint;
  typedef std::unordered_map<OnePoint2d,int,OnePoint2dHash>::iterator MyPointIt;
  typedef std::unordered_map<TwoPoint3d,int,TwoPoint3dHash> MyHash2Point;
  typedef std::unordered_map<TwoPoint3d,int,TwoPoint3dHash>::iterator My2PointIt;
#else
  typedef std::tr1::unordered_map<surfint,int> MySurfHash;
  typedef std::tr1::unordered_map<OnePoint2d,int,OnePoint2dHash> MyHashPoint;
  typedef std::tr1::unordered_map<OnePoint2d,int,OnePoint2dHash>::
    iterator MyPointIt;
  typedef std::tr1::unordered_map<TwoPoint3d,int,TwoPoint3dHash> MyHash2Point;
  typedef std::tr1::unordered_map<TwoPoint3d,int,TwoPoint3dHash>::
    iterator My2PointIt;
#endif

  MySurfHash *hash;           // hash for nlocal surf IDs
  int hashfilled;             // 1 if hash is filled with surf IDs

  Surf(class SPARTA *);
  ~Surf();
  void global(char *);
  void modify_params(int, char **);
  void init();
  void remove_ghosts();
  void add_line(surfint, int, double *, double *);
  void add_line_copy(int, Line *);
  void add_line_own(surfint, int, double *, double *);
  void add_line_temporary(surfint, int, double *, double *);
  void add_tri(surfint, int, double *, double *, double *);
  void add_tri_copy(int, Tri *);
  void add_tri_own(surfint, int, double *, double *, double *);
  void add_tri_own_clip(surfint, int, double *, double *, double *);
  void add_tri_temporary(surfint, int, double *, double *, double *);
  void rehash();
  void setup_owned();
  void setup_bbox();

  void compute_line_normal(int);
  void compute_tri_normal(int);
  void quad_corner_point(int, double *, double *, double *);
  void hex_corner_point(int, double *, double *, double *);

  double line_size(int);
  double line_size(Line *);
  double line_size(double *, double *);
  double axi_line_size(int);
  double axi_line_size(Line *);
  double tri_size(int, double &);
  double tri_size(Tri *, double &);
  double tri_size(double *, double *, double *, double &);

  void check_watertight_2d();
  void check_watertight_3d();
  void check_point_inside(int);
  void check_point_near_surf_2d();
  void check_point_near_surf_3d();

  void output_extent(int);
  double shortest_line(int);
  void smallest_tri(int, double &, double &);

  void add_collide(int, char **);
  int find_collide(const char *);
  void add_react(int, char **);
  int find_react(const char *);

  void group(int, char **);
  int add_group(const char *);
  int find_group(const char *);
  
  void compress_rebalance_explicit();
  void compress_rebalance_implicit();
  void reset_csurfs_implicit();

  void collate_vector(int, surfint *, double *, int, double *);
  void collate_vector_reduce(int, surfint *, double *, int, double *);
  void collate_vector_rendezvous(int, surfint *, double *, int, double *);

  void collate_array(int, int, surfint *, double **, double **);
  void collate_array_reduce(int, int, surfint *, double **, double **);
  void collate_array_rendezvous(int, int, surfint *, double **, double **);

  void redistribute_lines_clip(int, int);
  void redistribute_lines_temporary(int);
  void redistribute_tris_clip(int, int);
  void redistribute_tris_temporary(int);

  void write_restart(FILE *);
  void read_restart(FILE *);
  virtual void grow(int);
  virtual void grow_own(int);
  virtual void grow_temporary(int);
  bigint memory_usage();

 protected:
  int me,nprocs;
  int maxsc;                // max # of models in sc
  int maxsr;                // max # of models in sr
  
  // collate rendezvous data

  struct InRvousVec {
    surfint id;             // surface ID
    double value;           // compute value
  };

  double *out_rvous;
  int ncol_rvous;

  // watertight rendezvous data

  struct InRvousPoint {
    double x[2];            // 2d point coords
    int which;              // 1 for first endpoint, 2 for second endpoint
  };

  struct InRvousEdge {
    double x1[3],x2[3];     // 3d edge point coords
    int which;              // 1 for forward order, 2 for reverse order
  };

  // union data struct for packing 32-bit and 64-bit ints into double bufs
  // this avoids aliasing issues by having 2 pointers (double,int)
  //   to same buf memory
  // constructor for 32-bit int prevents compiler
  //   from possibly calling the double constructor when passed an int
  // copy to a double *buf:
  //   buf[m++] = ubuf(foo).d, where foo is a 32-bit or 64-bit int
  // copy from a double *buf:
  //   foo = (int) ubuf(buf[m++]).i;, where (int) or (tagint) match foo
  //   the cast prevents compiler warnings about possible truncation

  union ubuf {
    double d;
    int64_t i;
    ubuf(double arg) : d(arg) {}
    ubuf(int64_t arg) : i(arg) {}
    ubuf(int arg) : i(arg) {}
  };

  // private methods
  
  void point_line_compare(double *, double *, double *, double, int &, int &);
  void point_tri_compare(double *, double *, double *, double *, double *,
                         double, int &, int &, int, int, int);

  void collate_vector_allreduce(int, int *, double *, int, double *);
  void collate_vector_irregular(int, int *, double *, int, double *);
  void collate_array_allreduce(int, int, int *, double **, double **);
  void collate_array_irregular(int, int, int *, double **, double **);

  void check_watertight_2d_all();
  void check_watertight_2d_distributed();
  void check_watertight_3d_all();
  void check_watertight_3d_distributed();

  // callback functions for rendezvous communication
 
  static int rendezvous_vector(int, char *, int &, int *&, char *&, void *);
  static int rendezvous_array(int, char *, int &, int *&, char *&, void *);
  static int rendezvous_watertight_2d(int, char *, 
                                      int &, int *&, char *&, void *);
  static int rendezvous_watertight_3d(int, char *, 
                                      int &, int *&, char *&, void *);
  static int rendezvous_lines(int, char *, 
                              int &, int *&, char *&, void *);
  static int rendezvous_tris(int, char *, 
                             int &, int *&, char *&, void *);
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
