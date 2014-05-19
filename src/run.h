/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Cop2right (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(run,Run)

#else

#ifndef SPARTA_RUN_H
#define SPARTA_RUN_H

#include "pointers.h"

namespace SPARTA_NS {

class Run : protected Pointers {
 public:
  Run(class SPARTA *);
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Run command before grid is defined

UNDOCUMENTED

E: Run command before grid ghost cells are defined

UNDOCUMENTED

E: Invalid run command N value

The number of timesteps must fit in a 32-bit integer.  If you want to
run for more steps than this, perform multiple shorter runs.

E: Too many timesteps

The cummulative timesteps must fit in a 64-bit integer.

*/
