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

#ifdef FIX_CLASS

FixStyle(ablate,FixAblate)

#else

#ifndef SPARTA_FIX_ABLATE_H
#define SPARTA_FIX_ABLATE_H

#include "fix.h"

namespace SPARTA_NS {

class FixAblate : public Fix {
 public:
  int igroup;

  FixAblate(class SPARTA *, int, char **);
  ~FixAblate();
  int setmask();
  void init();
  void setup() {}
  void end_of_step();

  int pack_grid_one(int, char *, int);
  int unpack_grid_one(int, char *);
  void compress_grid();
  double memory_usage();

  void store_corners(int **);

 protected:
  int me;
  int groupbit,which,argindex,icompute,ifix,ncols;
  char *idsource;
  int storeflag;
   
  int nglocal;               // # of owned grid cells
  int nglocalmax;            // max size of per-cell vectors/arrays

  virtual void grow_percell(int);
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
