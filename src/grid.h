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

#ifdef SPARTA_MAP
#include <map>
#elif SPARTA_UNORDERED_MAP
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

#include "my_page.h"

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
  int maxsurfpercell;   // max surf elements in one child cell
  int maxlevel;         // max level of any child cell in grid, 0 = root
  int uniform;          // 1 if all child cells are at same level, else 0
  int unx,uny,unz;      // if uniform, effective global Nx,Ny,Nz of finest grid
  double cutoff;        // cutoff for ghost cells, -1.0 = infinite
  double cell_epsilon;  // half of smallest cellside of any cell in any dim
  int cellweightflag;   // 0/1+ for no/yes usage of cellwise fnum weighting

  int ngroup;               // # of defined groups
  char **gnames;            // name of each group
  int *bitmask;             // one-bit mask for each group

  // hash for all cell IDs (owned,ghost,parent)

#ifdef SPARTA_MAP
  std::map<cellint,int> *hash;
#elif SPARTA_UNORDERED_MAP
  std::unordered_map<cellint,int> *hash;
#else
  std::tr1::unordered_map<cellint,int> *hash;
#endif

  int hashfilled;             // 1 if hash is filled with parent/child IDs

  // list data structs

  MyPage<int> *csurfs;        // lists of surf indices for
                              // owned + ghost child cells with surfs
  MyPage<int> *csplits;       // lists of sub cell offsets for
                              // owned + ghost split info
  MyPage<int> *csubs;         // lists of sub cell indices for
                              // owned + ghost split info

  // owned or ghost child cell
  // includes unsplit cells, split cells, sub cells in any order
  // ghost cells are appended to owned

  struct ChildCell {          
    cellint id;               // cell ID in bitwise format
    int iparent;              // index of parent in pcells
    int proc;                 // proc that owns this cell
    int ilocal;               // index of this cell on owning proc
                              // must be correct for all kinds of ghost cells

    cellint neigh[6];         // info on 6 neighbor cells in cells/pcells
                              //   that fully overlap face
                              // order = XLO,XHI,YLO,YHI,ZLO,ZHI
                              // nmask stores flags for what all 6 represent
                              // if an index, store index into cells or pcells
                              // if unknown, store ID of neighbor child cell
                              // if non-periodic global boundary, ignored
    int nmask;                // 3 bits for each of 6 values in neigh
                              // 0 = index of child neighbor
                              // 1 = index of parent neighbor
                              // 2 = unknown child neighbor
                              // 3 = index of PBC child neighbor
                              // 4 = index of PBC parent neighbor
                              // 5 = unknown PBC child neighbor
                              // 6 = non-PBC boundary or ZLO/ZHI in 2d

    double lo[3],hi[3];       // opposite corner pts of cell
    int nsurf;                // # of surf elements in cell
                              // -1 = empty ghost cell
    int *csurfs;              // indices of surf elements in cell
                              // for sub cells, lo/hi/nsurf/csurfs
                              //   are same as in split cell containing them

    int nsplit;               // 1, unsplit cell
                              // N > 1, split cell with N sub cells
                              // N <= 0, neg of sub cell index (0 to Nsplit-1)
    int isplit;               // index into sinfo
                              // set for split and sub cells, -1 if unsplit
  };

  // info specific to owned child cell
  // includes unsplit cells, split cells, sub cells in same order as cells

  struct ChildInfo {
    int count;                // # of particles in this cell, 0 if split cell
    int first;                // index of 1st particle in this cell, -1 if none

    int mask;                 // grid group mask
    int type;                 // OUTSIDE,INSIDE,OVERLAP,UNKNOWN
    int corner[8];            // corner flags, 4/8 in 2d/3d
                              // OUTSIDE,INSIDE,OVERLAP,UNKNOWN
                              // ordered x first, y next, z last
                              // for sub cells, type/corner
                              //   are same as in split cell containing them

    double volume;            // flow volume of cell or sub cell
                              // entire cell volume for split cell
    double weight;            // fnum weighting for this cell
  };

  // additional info for owned or ghost split cell
  // ghost split cell info is appended to owned split cell info

  struct SplitInfo {
    int icell;                // index of split cell in cells this belongs to
    int xsub;                 // which sub cell (0 to Nsplit-1) xsplit is in
    double xsplit[3];         // coords of point in split cell
    int *csplits;             // sub cell (0 to Nsplit-1) each Nsurf belongs to
    int *csubs;               // indices in cells of Nsplit sub cells
  };

  // parent cell
  // global list of parent cells is stored by all procs

  struct ParentCell {
    cellint id;               // cell ID in bitwise format, 0 = root
    int mask;                 // grid group mask
    int level;                // level in hierarchical grid, 0 = root
    int nbits;                // # of bits to encode my ID, also my siblings
    int newbits;              // # of additional bits to encode my children
    int iparent;              // index of parent, -1 if id=root
    int grandparent;          // 1 if this cell is a grandparent, 0 if not
    int nx,ny,nz;             // sub grid within cell
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
  int nparent;                // # of parent cells

  int maxlocal;               // size of cinfo

  ChildCell *cells;           // list of owned and ghost child cells
  ChildInfo *cinfo;           // extra info for nlocal owned cells
  SplitInfo *sinfo;           // extra info for owned and ghost split cells
  ParentCell *pcells;         // global list of parent cells

  // restart buffers, filled by read_restart

  int nlocal_restart;
  cellint *id_restart;
  int *nsplit_restart;

  // methods

  Grid(class SPARTA *);
  ~Grid();
  void remove();
  void init();
  void add_child_cell(cellint, int, double *, double *);
  void add_parent_cell(cellint, int, int, int, int, double *, double *);
  void add_split_cell(int);
  void add_sub_cell(int, int);
  void remove_ghosts();
  void setup_owned();
  void acquire_ghosts();
  void rehash();
  void find_neighbors();
  void unset_neighbors();
  void reset_neighbors();
  void set_inout();
  void check_uniform();
  void type_check(int flag=1);
  void weight(int, char **);
  
  void refine_cell(int, int, int, int, int, int *, 
                   class Cut2d *, class Cut3d *);
  void coarsen_cell(int, int, int *, int *, int *, class AdaptGrid *,
                    class Cut2d *, class Cut3d *);

  void group(int, char **);
  int add_group(const char *);
  int find_group(const char *);

  void grow_pcells(int);

  void write_restart(FILE *);
  void read_restart(FILE *);
  int size_restart();
  int pack_restart(char *);
  int unpack_restart(char *);

  bigint memory_usage();

  void debug();

  // grid_comm.cpp

  int pack_one(int, char *, int, int, int);
  int unpack_one(char *, int, int);
  int pack_one_adapt(char *, char *, int);
  int pack_particles(int, char *, int);
  int unpack_particles(char *, int);
  void unpack_particles_adapt(int, char *);
  void compress();

  // grid_surf.cpp

  void surf2grid(int, int outflag=1);
  void surf2grid_one(int, int, int, int, class Cut3d *, class Cut2d *);
  void clear_surf();
  void clear_surf_restart();
  void combine_split_cell_particles(int, int);
  void assign_split_cell_particles(int);
  void allocate_surf_arrays();
  int *csubs_request(int);

  // grid_id.cpp

  int id_find_child(int, double *);
  int id_find_parent(cellint, cellint &);
  cellint id_str2num(char *);
  void id_num2str(int, char *);
  void id_pc_split(char *, char *, char *);
  void id_child_lohi(int, cellint, double *, double *);
  int id_bits(int, int, int);
  cellint id_find_face(double *, int, int, double *, double *);
  int id_child_from_parent_corner(int, int);

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

 private:
  int me;
  int maxcell;             // size of cells
  int maxsplit;            // size of sinfo
  int maxparent;           // size of pcells
  int maxbits;             // max bits allowed in a cell ID

  int neighmask[6];        // bit-masks for each face in nmask
  int neighshift[6];       // bit-shifts for each face in nmask

  // connection between one of my child cell corner pts
  // and a corner point in another child cell owned by me or another proc

  struct Connect {
    int newvalue;        // new value on the corner pt
    int icorner;         // what corner in other cell
    int icell;           // what other cell
    int proc;            // what proc owns the hcell
  };

  struct Connect2 {
    int newvalue;        // new value for neighbor cell
    int icell;           // index of other cell on owning proc
  };

  // bounding box for a clump of grid cells

  struct Box {
    double lo[3],hi[3];    // opposite corners of extended bbox
    int proc;              // proc that owns it
  };

  // Particle class values used for packing/unpacking particles in grid comm

  int ncustom;
  int nbytes_particle,nbytes_custom,nbytes_total;

  // private methods

  void acquire_ghosts_all();
  void acquire_ghosts_near();

  void box_intersect(double *, double *, double *, double *, 
                     double *, double *);
  int box_overlap(double *, double *, double *, double *);
  int box_periodic(double *, double *, Box *);

  void grow_cells(int, int);
  void grow_sinfo(int);

  void surf2grid_stats();
  void flow_stats();
  double flow_volume();

  // grid_comm.cpp
  // static variable for ring communication callback to access class data
  // callback function for ring communication

  static Grid *gptr;
  static void unpack_ghosts(int, char *);
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
