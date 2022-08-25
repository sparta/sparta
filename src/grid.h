/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_GRID_H
#define SPARTA_GRID_H

#include "stdio.h"
#include "pointers.h"
#include "hash3.h"
#include "my_page.h"
#include "surf.h"

namespace SPARTA_NS {

class Grid : protected Pointers {
 public:
  int exist;            // 1 if grid is defined
  int exist_ghost;      // 1 if ghost cells exist
  int clumped;          // 1 if grid ownership is clumped, due to RCB
                        // if not, some operations are not allowed

  bigint ncell;         // global count of child cells (unsplit+split, no sub)
  bigint nunsplit;      // global count of unsplit cells
  int nsplit;           // global count of split cells
  int nsub;             // global count of split sub cells

  int maxlevel;         // max level of any child cell in grid, 0 = root
  int plevel_limit;     // allocation bound of plevels

  int uniform;          // 1 if all child cells are at same level, else 0
  int unx,uny,unz;      // if uniform, effective global Nx,Ny,Nz of finest grid
  double cutoff;        // cutoff for ghost cells, -1.0 = infinite
  double cell_epsilon;  // half of smallest cellside of any cell in any dim
  int cellweightflag;   // 0/1+ for no/yes usage of cellwise fnum weighting

  int surfgrid_algorithm;  // algorithm for overlap of surfs & grid cells
  int maxsurfpercell;   // max surf elements in one child cell
  int maxsplitpercell;  // max split cells in one child cell

  double dt_global;     // global timestep size
  double time_global;   // the global time
  bool use_cell_dt = false;

  int ngroup;           // # of defined groups
  char **gnames;        // name of each group
  int *bitmask;         // one-bit mask for each group
  int *inversemask;     // inverse mask for each group

  double tmap,tsplit;   // timing breakdowns of both grid2surf() algs
  double tcomm1,tcomm2,tcomm3,tcomm4;

  int copy,copymode;    // 1 if copy of class (prevents deallocation of
                        //  base class when child copy is destroyed)

  // cell ID hash (owned + ghost, no sub-cells)

#ifdef SPARTA_MAP
  typedef std::map<cellint,int> MyHash;
#elif SPARTA_UNORDERED_MAP
  typedef std::unordered_map<cellint,int> MyHash;
#else
  typedef std::tr1::unordered_map<cellint,int> MyHash;
#endif

  MyHash *hash;
  int hashfilled;             // 1 if hash is filled with cell IDs

  // list data structs

  MyPage<surfint> *csurfs;    // lists of surf indices for
                              // owned + ghost child cells with surfs
  MyPage<int> *csplits;       // lists of sub cell offsets for
                              // owned + ghost split info
  MyPage<int> *csubs;         // lists of sub cell indices for
                              // owned + ghost split info

  // owned or ghost child cell
  // includes unsplit cells, split cells, sub cells in any order
  // ghost cells are appended to owned

  struct SPARTA_ALIGN(128) ChildCell {
    cellint id;               // ID of child cell
    int level;                // level of cell in hierarchical grid, 0 = root
    int proc;                 // proc that owns this cell
    int ilocal;               // index of this cell on owning proc
                              // must be correct for all ghost cells

    cellint neigh[6];         // info on 6 neighbor cells that fully overlap faces
                              // order = XLO,XHI,YLO,YHI,ZLO,ZHI
                              // nmask has flags for what the values represent
                              // for nmask 0,3: index into cells (own or ghost)
                              // for nmask 1,4: index into pcells
                              // for nmask 2,5: neigh is unknown, value ignored
                              // if unknown or non-PBC boundary or 2d ZLO/ZHI, ignored
                              // must be cellint, b/c sometimes stores cell IDs

    int nmask;                // 3 bits for each of 6 values in neigh
                              // 0 = index of child neighbor I own
                              // 1 = index of parent neighbor in pcells
                              // 2 = unknown child neighbor
                              // 3 = index of PBC child neighbor I own
                              // 4 = index of PBC parent neighbor in pcells
                              // 5 = unknown PBC child neighbor
                              // 6 = non-PBC boundary or ZLO/ZHI in 2d

    double lo[3],hi[3];       // opposite corner pts of cell

    int nsurf;                // # of surf elements in cell
                              // -1 = empty ghost cell
    surfint *csurfs;          // indices of surf elements in cell
                              // sometimes global surf IDs are stored
                              // for sub cells, lo/hi/nsurf/csurfs
                              //   are same as in split cell containing them

    int nsplit;               // 1, unsplit cell
                              // N > 1, split cell with N sub cells
                              // N <= 0, neg of sub cell index (0 to Nsplit-1)
    int isplit;               // index into sinfo
                              // set for split and sub cells, -1 if unsplit

    double dt_desired;        // desired timestep for this cell
    double time;              // time for this cell
  };

  // info specific to owned child cell
  // includes unsplit cells, split cells, sub cells in same order as cells

  struct ChildInfo {
    int count;                // # of particles in this cell, 0 if split cell
    int first;                // index of 1st particle in this cell, -1 if none

    int mask;                 // grid group mask
    int type;                 // OUTSIDE,INSIDE,OVERLAP,UNKNOWN
    int corner[8];            // corner flags, 4/8 in 2d/3d
                              // OUTSIDE,INSIDE,UNKNOWN
                              // no OVERLAP is used for this anymore I think
                              // ordered x first, y next, z last
                              // for sub cells, type/corner
                              //   are same as in split cell containing them

    double volume;            // flow volume of cell or sub cell
                              // entire cell volume for split cell
    double weight;            // fnum weighting for this cell
  };

  // additional info for owned or ghost split cell or sub cell
  // ghost split cell info is appended to owned split cell info

  struct SplitInfo {
    int icell;                // index of split cell this sub cell belongs to
    int xsub;                 // which sub cell (0 to Nsplit-1) xsplit is in
    double xsplit[3];         // coords of point in split cell
    int *csplits;             // sub cell (0 to Nsplit-1) each Nsurf belongs to
    int *csubs;               // indices in cells of Nsplit sub cells
  };

  struct ParentLevel {
    int nbits;                // nbits = # of bits to store parent ID at this level
    int newbits;              // newbits = extra bits to store children of this parent
    int nx,ny,nz;             // child grid below this parent level
    bigint nxyz;              // # of child cells of this parent level
  };

  struct ParentCell {
    cellint id;               // ID of parent cell
    double lo[3],hi[3];       // opposite corner pts of cell
  };

  int nlocal;                 // # of child cells I own (all 3 kinds)
  int nghost;                 // # of ghost child cells I store (all 3 kinds)
  int nempty;                 // # of empty ghost cells I store
  int nunsplitlocal;          // # of unsplit cells I own
  int nunsplitghost;          // # of ghost unsplit cells I store
  int nsplitlocal;            // # of split cells I own
  int nsplitghost;            // # of ghost split cells I store
  int nsublocal;              // # of sub cells I own
  int nsubghost;              // # of ghost sub cells I store
  int nparent;                // # of parent cell neighbors I store

  int maxlocal;               // size of cinfo
  int maxparent;              // size of pcells

  ChildCell *cells;           // list of owned and ghost child cells
  ChildInfo *cinfo;           // extra info for nlocal owned cells
  SplitInfo *sinfo;           // extra info for owned and ghost split cells

  ParentLevel *plevels;       // list of parent levels, level = root = simulation box
  ParentCell *pcells;         // list of parent cell neighbors

  // restart buffers, filled by read_restart

  int nlocal_restart;
  cellint *id_restart;
  int *level_restart,*nsplit_restart;
  double *dt_desired_restart;
  double *time_restart;

  // methods

  Grid(class SPARTA *);
  ~Grid();
  void remove();
  void init();
  void add_child_cell(cellint, int, double *, double *);
  void add_split_cell(int);
  void add_sub_cell(int, int);
  void notify_changed();
  int set_minlevel();
  void set_maxlevel();
  void setup_owned();
  void remove_ghosts();
  void acquire_ghosts(int surfflag=1);
  void rehash();
  void find_neighbors();
  void unset_neighbors();
  void reset_neighbors();
  void set_inout();
  void check_uniform();
  void type_check(int flag=1);
  void weight(int, char **);
  void weight_one(int);

  void refine_cell(int, int *, class Cut2d *, class Cut3d *);
  void coarsen_cell(cellint, int, double *, double *,
                    int, int *, int *, int *, void **, char **,
                    class Cut2d *, class Cut3d *);

  void group(int, char **);
  int add_group(const char *);
  int find_group(const char *);
  int check_uniform_group(int, int *, double *, double *);

  void write_restart(FILE *);
  void read_restart(FILE *);
  int size_restart();
  int size_restart(int);
  int pack_restart(char *);
  int unpack_restart(char *);

  bigint memory_usage();

  void debug();

  // grid_comm.cpp

  int pack_one(int, char *, int, int, int, int);
  int unpack_one(char *, int, int, int, int sortflag=0);
  int pack_one_adapt(char *, char *, int);
  int pack_particles(int, char *, int);
  int unpack_particles(char *, int, int);
  void unpack_particles_adapt(int, char *);
  void compress();

  // grid_surf.cpp

  void surf2grid(int, int outflag=1);
  void surf2grid_implicit(int, int outflag=1);
  void surf2grid_one(int, int, int, int, class Cut3d *, class Cut2d *);
  void clear_surf();
  void clear_surf_restart();
  void combine_split_cell_particles(int, int);
  void assign_split_cell_particles(int);
  int outside_surfs(int, double *, class Cut3d *, class Cut2d *);
  void allocate_surf_arrays();
  int *csubs_request(int);

  // grid_id.cpp

  void id_point_child(double *, double *, double *, int, int, int,
                      int &, int &, int &);
  cellint id_parent_of_child(cellint, int);
  int id_find_child(cellint, int, double *, double *, double *);
  cellint id_uniform_level(int, int, int, int);
  void id_find_child_uniform_level(int, int, double *, double *, double *,
                                   int &, int &, int &);
  cellint id_neigh_same_parent(cellint, int, int);
  cellint id_neigh_same_level(cellint, int, int);
  cellint id_refine(cellint, int, int);
  cellint id_coarsen(cellint, int);
  cellint id_ichild(cellint, cellint, int);
  int id_level(cellint);
  void id_child_lohi(int, double *, double *, cellint, double *, double *);
  void id_lohi(cellint, int, double *, double *, double *, double *);
  int id_bits(int, int, int);
  void id_num2str(cellint, char *);

  // extract/return neighbor flag for iface from per-cell nmask
  // inlined for efficiency

  inline int neigh_decode(int nmask, int iface) {
    return (nmask & neighmask[iface]) >> neighshift[iface];
  }

  // overwrite neighbor flag for iface in per-cell nmask
  // first line zeroes the iface bits via one's complement of mask
  // inlined for efficiency
  // return updated nmask

  inline int neigh_encode(int flag, int nmask, int iface) {
    nmask &= ~neighmask[iface];
    nmask |= flag << neighshift[iface];
    return nmask;
  }

  double get_particle_time(cellint, double); // get particle time as a function of cell desired dt, global_time, and RNG


 protected:
  int me;
  int maxcell;             // size of cells
  int maxsplit;            // size of sinfo
  int maxbits;             // max bits allowed in a cell ID

  int neighmask[6];        // bit-masks for each face in nmask
  int neighshift[6];       // bit-shifts for each face in nmask

  class Cut2d *cut2d;
  class Cut3d *cut3d;

  // connection between one of my cells and a neighbor cell on another proc

  struct Connect {
    int itype;           // type of sending cell
    int marktype;        // new type value (IN/OUT) for neighbor cell
    int jcell;           // index of neighbor cell on receiving proc (owner)
  };

  // bounding box for a clump of grid cells

  struct Box {
    double lo[3],hi[3];    // opposite corners of extended bbox
    int proc;              // proc that owns it
  };

  // tree of RCB cuts for a partitioned uniform block of grid cells

  struct GridTree {
    int dim,cut;
  };

  // data structs for rendezvous comm

  struct InRvous {
    int proc;
    surfint surfID;
  };

  struct OutRvousLine {
    Surf::Line line;
  };

  struct OutRvousTri {
    Surf::Tri tri;
  };

  // surf ID hashes

#ifdef SPARTA_MAP
    typedef std::map<surfint,int> MySurfHash;
    typedef std::map<surfint,int>::iterator MyIterator;
#elif SPARTA_UNORDERED_MAP
    typedef std::unordered_map<surfint,int> MySurfHash;
    typedef std::unordered_map<surfint,int>::iterator MyIterator;
#else
    typedef std::tr1::unordered_map<surfint,int> MySurfHash;
    typedef std::tr1::unordered_map<surfint,int>::iterator MyIterator;
#endif

  // Particle class values used for packing/unpacking particles in grid comm

  int ncustom;
  int nbytes_particle,nbytes_custom,nbytes_total;

  // private methods

  void surf2grid_cell_algorithm(int);
  void surf2grid_surf_algorithm(int);
  void surf2grid_split(int, int);
  void recurse2d(cellint, int, double *, double *,
                 int, Surf::Line *, double *, double *,
                 int &, int &, int **&, MyHash *, MyHash *);
  void recurse3d(cellint, int, double *, double *,
                 int, Surf::Tri *, double *, double *,
                 int &, int &, int **&, MyHash *, MyHash *);
  void partition_grid(int, int, int, int, int, int, int, int, GridTree *);
  void mybox(int, int, int, int &, int &, int &, int &, int &, int &,
             GridTree *);
  void box_drop(int *, int *, int, int, GridTree *, int &, int *);

  void acquire_ghosts_all(int);
  void acquire_ghosts_near(int);
  void acquire_ghosts_near_less_memory(int);

  void box_intersect(double *, double *, double *, double *,
                     double *, double *);
  int box_overlap(double *, double *, double *, double *);
  int box_periodic(double *, double *, Box *);

  virtual void grow_cells(int, int);
  virtual void grow_pcells();
  virtual void grow_sinfo(int);

  void surf2grid_stats();
  void flow_stats();
  double flow_volume();

  // callback function for ring communication and class variable to access

  int unpack_ghosts_surfflag;
  static void unpack_ghosts(int, char *, void *);

  // callback functions for rendezvous communication

  static int rendezvous_surfrequest(int, char *, int &, int *&, char *&, void *);
};

}

#endif

/* ERROR/WARNING messages:

E: Cell ID has too many bits

Cell IDs must fit in 32 bits (SPARTA small integer) or 64 bits (SPARTA
big integer), as specified by the -DSPARTA_SMALL, -DSPARTA_BIG, or
-DSPARTA_BIGBIG options in the machine Makefile used to build
SPARTA.  See Section 2.2 of the manual for details.  And see Section
4.8 for details on how cell IDs are formatted.

E: Owned cells with unknown neighbors = %d

One or more grid cells have unknown neighbors which will prevent
particles from moving correctly.  Please report the issue to the
SPARTA developers.

E: Grid in/out self-mark error %d for icell %d, icorner %d, connect %d %d, other cell %d, other corner %d, values %d %d\n

A grid cell was incorrectly marked as inside, outside, or overlapping
with surface elements.  Please report the issue to the SPARTA
developers.

E: Grid in/out other-mark error %d\n

Grid cell marking as inside, outside, or overlapping with surface
elements failed.  Please report the issue to the SPARTA developers.

E: Cell type mis-match when marking on self

Grid cell marking as inside, outside, or overlapping with surface
elements failed.  Please report the issue to the SPARTA developers.

E: Parent cell child missing

Hierarchical grid traversal failed.  Please report the issue to the
SPARTA developers.

E: Cell type mis-match when marking on neigh proc

Grid cell marking as inside, outside, or overlapping with surface
elements failed.  Please report the issue to the SPARTA developers.

E: Grid cells marked as unknown = %d

Grid cell marking as inside, outside, or overlapping with surface
elements did not successfully mark all cells.  Please report the issue
to the SPARTA developers.

W: Grid cell interior corner points marked as unknown = %d

Corner points of grid cells interior to the simulation domain were not
all marked successfully as inside, outside, or overlapping with
surface elements.  This should normally not happen, but does
not affect simulations.

E: Grid cell corner points on boundary marked as unknown = %d

Corner points of grid cells on the boundary of the simulation domain
were not all marked successfully as inside, outside, or overlapping
with surface elements.  Please report the issue to the SPARTA
developers.

E: Cannot weight cells before grid is defined

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Cannot use weight cell radius unless axisymmetric

An axisymmetric model is required for this style of cell weighting.

*/
