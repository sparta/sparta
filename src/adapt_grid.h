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

#ifdef COMMAND_CLASS

CommandStyle(adapt_grid,AdaptGrid)

#else

#ifndef SPARTA_ADAPT_GRID_H
#define SPARTA_ADAPT_GRID_H

#include "pointers.h"
#include "hash3.h"

namespace SPARTA_NS {

class AdaptGrid : protected Pointers {
 public:
  int mode;                 // 0 = adapt_grid command, 1 = fix adapt command
  int action1,action2;      // adaptation actions, so FixAdapt can access
  char *file;               // output file name, so FixAdapt can access

  struct SendAdapt {
    cellint parentID;       // parentID this child cell info is for
    int plevel;             // grid level of parentID
    int owner;              // proc which owns the parent cell for coarsening
    int proc;               // proc which owns this child cell
    int icell;              // local index of child cell
    int type;               // type of child cell
    int ichild;             // which child within parent cell (0 to Nxyz-1)
    int nsurf;              // # of surfs in child cell
    int np;                 // # of particles in child cell or all its sub cells
  };

  AdaptGrid(class SPARTA *);
  ~AdaptGrid();
  void command(int, char **);
  void process_args(int, char **);
  void check_args(int);
  void setup(int);
  bigint refine();
  bigint coarsen();
  void cleanup();
  void write_file();

 private:
  int me,nprocs;
  int groupbit,style,niterate,sgroupbit,regstyle;
  int rdecide,cdecide,combine,minlevel,maxlevel;
  int nx,ny,nz;
  int valuewhich,valindex,icompute,ifix;
  double rcount,ccount,rvalue,cvalue,rfrac,cfrac;
  double surfsize;
  double sdir[3];
  char *computeID,*valueID;
  class Region *region;
  class Compute *compute;
  class Fix *fix;

  int *childlist;
  class RanKnuth *random;
  class Cut3d *cut3d;
  class Cut2d *cut2d;

  // rlist = list of child cells I will refine

  int *rlist;
  int rnum;

  // clist = list of parent cells I may coarsen
  // action = list of parent cells I will coarsen

  struct CList {
    cellint parentID;       // ID of parent cell
    int flag;               // 1 if coarsen, 0 if not
    int plevel;             // level of parent cell
    int nchild;             // # of children of parent cell
    int nexist;             // # of identified children (owned or recvd)
    int *proc;              // proc that owns each child
    int *index;             // local index of each child on owning proc
    double *value;          // style value for each child
  };

  struct ActionList {
    cellint parentID;       // ID of parent cell
    int plevel;             // level of parent cell
    int anyinside;          // 1 if any children are an INSIDE cell
    int nchild;             // # of children of parent cell
    int *index;             // local icell if I own each child cell, else -1
    int *nsurf;             // # of surfs in each child cell
    int *np;                // # of particles in each child cell
    void **surfs;           // list of surf info per cell (indices or lines/tris)
    char **particles;       // list of particles in each child cell
  };

  struct CList *clist;
  struct ActionList *alist;
  int cnum,cnummax;
  int anum,anummax;

  char *spbuf;                // buf allocated by Comm::send_cells_adapt
                              // for child cell surf and particle data

  // hash for child cell IDs created by refinement or coarsening

#ifdef SPARTA_MAP
  typedef std::map<cellint,int> MyHash;
#elif SPARTA_UNORDERED_MAP
  typedef std::unordered_map<cellint,int> MyHash;
#else
  typedef std::tr1::unordered_map<cellint,int> MyHash;
#endif

  MyHash *chash,*rhash;

  // Rvous1 send info

  struct Rvous1 {
    cellint parentID;
    int plevel,proc,icell,ichild;     // ichild = 0 to Nxyz-1
    double value;
  };

  // Rvous2 send info

  struct Rvous2 {
    cellint parentID;
    int owner,icell,ichild;
  };

  // methods

  void candidates_refine();
  void refine_particle();
  void refine_surf();
  void refine_value();
  void refine_random();
  int perform_refine();

  void candidates_coarsen();
  double coarsen_particle_cell(int);
  double coarsen_surf_cell(int);
  double coarsen_value_cell(int);
  void coarsen_particle();
  void coarsen_surf();
  void coarsen_value();
  void coarsen_random();
  void particle_surf_comm();
  int perform_coarsen();

  double value_compute(int);
  double value_fix(int);

  int region_check(double *, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
