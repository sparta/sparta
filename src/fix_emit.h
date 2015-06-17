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

#ifndef SPARTA_FIX_EMIT_H
#define SPARTA_FIX_EMIT_H

#include "fix.h"

namespace SPARTA_NS {

class FixEmit : public Fix {
 public:
  FixEmit(class SPARTA *, int, char **);
  virtual ~FixEmit();
  int setmask();
  virtual void init();
  void start_of_step();
  double compute_vector(int);

  void add_grid_one(int, int);
  int pack_grid_one(int, char *, int);
  int unpack_grid_one(int, char *);
  void compress_grid();
  virtual void post_compress_grid() {}

 protected:
  int perspecies;
  class Region *region;
  class RanPark *random;
  int nsingle,ntotal;

  int nglocal;         // copy of cell->nlocal
  int nglocalmax;      // max size of c2list
  int *c2list;         // index into clist for each owned cell
                       // -1 if no tasks for the cell
                       // only for unsplit and split cells

  int nlist;           // # of owned cells with insert tasks
  int nlistmax;        // max size of clist,clistnum,clistfirst
  int *clist;          // local indices of cells with insert tasks
  int *clistnum;       // # of insert tasks in each cell with tasks
  int *clistfirst;     // first task ID for each cell with tasks

  int ntask;           // # of insert tasks in underlying child class

  virtual int create_task(int) = 0;
  virtual void perform_task() = 0;
  virtual int pack_task(int, char *, int) = 0;
  virtual int unpack_task(char *, int) = 0;
  virtual void copy_task(int, int, int, int) = 0;

  void grow_percell(int);
  void grow_list();
  double mol_inflow(double, double, double);
  void options(int, char **);
  virtual int option(int, char **);
};

}

#endif

/* ERROR/WARNING messages:

*/
