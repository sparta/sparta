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

FixStyle(ave/surf,FixAveSurf)

#else

#ifndef SPARTA_FIX_AVE_SURF_H
#define SPARTA_FIX_AVE_SURF_H

#include "fix.h"
#include "hash3.h"

namespace SPARTA_NS {

class FixAveSurf : public Fix {
 public:
  FixAveSurf(class SPARTA *, int, char **);
  ~FixAveSurf();
  int setmask();
  void init();
  void setup();
  void end_of_step();
  double memory_usage();

 private:
  int groupbit;
  int nvalues,maxvalues;
  int nrepeat,irepeat,nsample,ave;
  bigint nvalid;
  int *which,*argindex,*value2index;
  char **ids;

  int nsurf;               // # of explicit surfs I store = surf->nlocal+nghost
  int nown;                // # of explicit surfs I own = surf->nown
  double *accvec;          // accumulation vector
  double **accarray;       // accumulation array
  double *bufvec;          // surf collate vector for surfs I own
  double **bufarray;       // surf collate array for surfs I own
  int *masks;              // surface element group masks for surfs I own

  int ntally;              // # of surfs I have tallies for
  int maxtally;            // # of tallies currently allocated
  surfint *tally2surf;     // tally2surf[I] = surf ID of Ith tally
  double *vec_tally;       // tally values, maxtally in length
  double **array_tally;

  // hash for surf IDs

#ifdef SPARTA_MAP
  typedef std::map<surfint,int> MyHash;
#elif defined SPARTA_UNORDERED_MAP
  typedef std::unordered_map<surfint,int> MyHash;
#else
  typedef std::tr1::unordered_map<surfint,int> MyHash;
#endif

  MyHash *hash;

  void options(int, int, char **);
  void grow_tally();
  bigint nextvalid();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute ID for fix ave/surf does not exist

Self-explanatory.

E: Fix ave/surf compute does not calculate per-surf values

Self-explanatory.

E: Fix ave/surf compute does not calculate a per-surf vector

Self-explanatory.

E: Fix ave/surf compute does not calculate a per-surf array

Self-explanatory.

E: Fix ave/surf compute array is accessed out-of-range

Self-explanatory.

E: Fix ID for fix ave/surf does not exist

Self-explanatory.

E: Fix ave/surf fix does not calculate per-surf values

Self-explanatory.

E: Fix ave/surf fix does not calculate a per-surf vector

Self-explanatory.

E: Fix ave/surf fix does not calculate a per-surf array

Self-explanatory.

E: Fix ave/surf fix array is accessed out-of-range

Self-explanatory.

E: Fix for fix ave/surf not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/surf is
requesting a value on a non-allowed timestep.

E: Variable name for fix ave/surf does not exist

Self-explanatory.

E: Fix ave/surf variable is not surf-style variable

Self-explanatory.

*/
