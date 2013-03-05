/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_CUT2D_H
#define SPARTA_CUT2D_H

#include "pointers.h"
#include "my_vec.h"
#include "my_double_linked_list.h"
#include "my_page.h"

namespace SPARTA_NS {

class Cut2d : protected Pointers {
 public:
  Cut2d(class SPARTA *);
  ~Cut2d();
  void surf2grid();
  void split();

 private:
  // just for VERBOSE output
  int icell;

  MyVec<int> used;
  MyVec<int> startpts;
  MyVec<int> endpts;
  
  struct PLone {
    int iline;
    PLone *prev,*next;
  };

  struct Cpt {
    int flag;
    double x[2];
    double dot;
    int ipl,oindex;
    Cpt *prev,*next;
  };

  struct Opt {
    int flag;
    double x[2];
    int cindex;
  };

  struct Entrypt {
    int iopt;
    int index;
    Entrypt *prev,*next;
  };
  
  struct Loop {
    int flag;
    double area;
    MyVec<int> lines;
  };

  struct PG {
    double area;
    MyVec<int> lines;
  };

  MyPage<PLone> ppool;
  MyVec< MyDoubleLinkedList<PLone*> > pl;

  MyVec< MyVec<Opt> > opts;

  MyPage<Cpt> cpool;
  MyDoubleLinkedList<Cpt*> cpts;
  Cpt *cindex[5];

  MyPage<Entrypt> epool;
  MyDoubleLinkedList<Entrypt*> entrypts;

  MyVec<Loop> loops;
  MyVec<PG> pg;

  void line2pl(int, int *);
  void weiler_intersect(double *, double *, int *);
  void interleave(double *, double *, double *, double *, int, int, int);
  void weiler_walk(int, double *, double *, int *);
  void loop2pg(int, double *, double *, int *);
  void surf2pg(int, int *, int *);
  int split_point(double *, double *, int, int *, int *, double *);

  int cliptest(double *, double *, double *, double *);
  void clip(double *, double *, double *, double *, double *, double *);
  int sameborder(double *pt1, double *pt2, double *lo, double *hi);
  int corner(double *, double *, double *);
  int ptflag(double *pt, double *lo, double *hi);
};

}

#endif

/* ERROR/WARNING messages:

E: Bad grid of processors for create_grid

UNDOCUMENTED

E: Per-processor grid count is too big

UNDOCUMENTED

*/
