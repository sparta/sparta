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

#ifdef SPARTA_MAP
#include <map>
#elif SPARTA_UNORDERED_MAP
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

namespace SPARTA_NS {

class AdaptGrid : protected Pointers {
 public:
  int action1,action2;      // adaptation actions, so FixAdapt can access

  struct SendAdapt {
    int proc;               // proc that is sending
    int icell;              // index of my cell to send
    int ichild;             // which child cell (0 to NxNyNz-1) within parent
    int iparent;            // index of cell's parent cell
    int type;               // cinfo type of cell
    int nsurf;              // # of surfs in cell
    int np;                 // # of particles in cell or all its sub cells
    double value;           // per-cell value to coarsen on
  };

  AdaptGrid(class SPARTA *);
  ~AdaptGrid();
  void command(int, char **);
  void process_args(int, int, char **);
  void check_args(int, int);
  void setup();
  int refine();
  int coarsen(int);
  void cleanup();

 private:
  int me,nprocs;
  int style,niterate;
  char *computeID;
  double rthresh,cthresh,cvalue,rfrac,cfrac;
  int minlevel,maxlevel;
  int nx,ny,nz;
  class Region *region;
  int rstyle;
  char *valueID;
  int valuewhich,valueindex;
  int icompute,ifix,argindex;

  int *newcell;
  class RanPark *random;
  int *delparent;
  class Cut3d *cut3d;
  class Cut2d *cut2d;

  // rlist = list of child cells I will refine

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

  // sadapt = list of child cell info to send to other procs for coarsening

  struct SendAdapt *sadapt;
  int nsend,nrecv,replyany;

  SendAdapt **sa_header;
  int **sa_csurfs;
  char **sa_particles;

  // hash for child cell IDs created by converting a coarsened parent cell

#ifdef SPARTA_MAP
  std::map<cellint,int> *chash;
#elif SPARTA_UNORDERED_MAP
  std::unordered_map<cellint,int> *chash;
#else
  std::tr1::unordered_map<cellint,int> *chash;
#endif

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
  void coarsen_particle();
  void coarsen_surf();
  void coarsen_value();
  void coarsen_random();
  int perform_coarsen();
  void gather_parents_coarsen(int, int);

  int region_check(double *, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
