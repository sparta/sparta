/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(lammps,FixLAMMPS)

#else

#ifndef SPARTA_FIX_LAMMPS_H
#define SPARTA_FIX_LAMMPS_H

#include "fix.h"
#include "grid.h"
#include "surf.h"
#include <string>

namespace SPARTA_NS {

class FixLAMMPS : public Fix {
 public:
  FixLAMMPS(class SPARTA *, int, char **);
  ~FixLAMMPS();
  int setmask();
  void init();
  void end_of_step();

 protected:
  int ncells, dimension;
  int groupbit;
  int numPar, fnum;
  double KE_ave, Tave, p_spa;
  double *v, mass, ***v_lmp;
  bigint nvalid;
  FILE *fp;                  // file to write dump to
  char filename[256];           // LAMMPS input file
  char filenameS[256];          // file specific to proc and cell
  bigint nextvalid();

  // Surface Temperature Specific
  int firstflag, tindex, icompute;
  double *tvector_me, *tvector;
  double twall;
  char *id_nrho;
  class Compute *c_nrho;

  class RanKnuth *random;     // RNG for particle reflection
  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Compute ID for fix lammps does not exist

Double-check ID used for fix.

E: 

*/