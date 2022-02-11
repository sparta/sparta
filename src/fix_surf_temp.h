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

FixStyle(surf/temp,FixSurfTemp)

#else

#ifndef SPARTA_FIX_SURF_TEMP_H
#define SPARTA_FIX_SURF_TEMP_H

#include "fix.h"
#include "surf.h"

namespace SPARTA_NS {

class FixSurfTemp : public Fix {
 public:
  FixSurfTemp(class SPARTA *, int, char **);
  ~FixSurfTemp();
  int setmask();
  void init();
  virtual void end_of_step();

 private:
  int ifix,m;
  double twall,prefactor,emi;
  char *id_qw;
  int qwindex;
  class Fix *fqw;
  double *qw,*qw_all;

  Surf::Line *lines;
  Surf::Tri *tris;

  int distributed,implicit,dimension;  // Surf settings
  int firstflag;

  int groupbit;              // mask for surface group
  int nown;                  // # of surf elements owned by this proc
  int nchoose;               // # of surf elements output by this proc
  int nsurf;                 // nown or nlocal
  int *cglobal;              // indices of global elements for nchoose
  int *clocal;               // indices of local owned elements for nchoose
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Cannot use non-rcb fix balance with a grid cutoff

This is because the load-balancing will generate a partitioning
of cells to processors that is dispersed and which will not work
with a grid cutoff >= 0.0.

*/
