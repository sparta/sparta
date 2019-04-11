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
    int proc;               // proc that is sending
    int icell;              // index of my cell to send
    int ichild;             // which child cell (0 to NxNyNz-1) within parent
    int iparent;            // index of cell's parent cell
    int type;               // cinfo type of cell
    int mask;               // cinfo mask of cell
    int nsurf;              // # of surfs in cell
    int np;                 // # of particles in cell or all its sub cells
    double value;           // per-cell value to coarsen on
  };

  SendAdapt **sa_header;    // communicated coarsening info, so Grid can access
  surfint **sa_csurfs;
  char **sa_particles;

  AdaptGrid(class SPARTA *);
  ~AdaptGrid();
  void command(int, char **);
  void process_args(int, char **);
  void check_args(int);
  void setup(int);
  int refine();
  int coarsen(int);
  void add_grid_fixes();
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
  class RanPark *random;
  int *delparent;
  class Cut3d *cut3d;
  class Cut2d *cut2d;

  // rlist = list of child cells I may refine

  int *rlist;
  int rnum;

  // parent cell info for coarsening

  int *pcount,*powner;

  // ctask = list of coarsen tasks I will perform on parent cells

  struct CTask {
    cellint id;             // id of parent cell
    int iparent;            // index of parent cell
    int flag;               // 1 if coarsen, 0 if not
    int nchild;             // # of children of parent cell
    int ncomplete;          // # which info is complete for (owned or recvd)
    int *proc;              // proc that owns each child
    int *index;             // local index of each child on owning proc
    int *recv;              // index into recv list of children from other procs
  };

  struct CTask *ctask;
  int cnum,cnummax;

  // newcells = IDs of new child cells I create, either via refine or coarsen

  cellint *newcells;
  int nnew,maxnew;

  // sadapt = list of child cell info to send to other procs for coarsening

  struct SendAdapt *sadapt;
  int nsend,nrecv,replyany;

  // hash for child cell IDs created by converting a coarsened parent cell

#ifdef SPARTA_MAP
  typedef std::map<cellint,int> MyHash;
#elif SPARTA_UNORDERED_MAP
  typedef std::unordered_map<cellint,int> MyHash;
#else
  typedef std::tr1::unordered_map<cellint,int> MyHash;
#endif

  MyHash *chash;

  // methods

  void candidates_refine();
  void refine_particle();
  void refine_surf();
  void refine_value();
  void refine_random();
  int perform_refine();
  void gather_parents_refine(int, int);

  void assign_parents_coarsen(int);
  void candidates_coarsen(int);
  void coarsen_group();
  void coarsen_particle();
  void coarsen_surf();
  void coarsen_value();
  void coarsen_random();
  int perform_coarsen();
  void gather_parents_coarsen(int, int);

  double value_compute(int);
  double value_fix(int);

  int region_check(double *, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
