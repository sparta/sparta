/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
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

Self-explanatory.

E: Run command before grid ghost cells are defined

Normally, ghost cells will be defined when the grid is created via the
create_grid or read_grid commands.  However, if the global gridcut
cutoff is set to a value >= 0.0, then ghost cells can only be defined
if the partiioning of cells to processors is clumped, not dispersed.
See the fix balance command for an explanation.  Invoking the fix
balance command with a clumped option will trigger ghost cells to be
defined.

E: Invalid run command N value

The number of timesteps must fit in a 32-bit integer.  If you want to
run for more steps than this, perform multiple shorter runs.

E: Invalid run command upto value

Self-explanatory.

E: Invalid run command start/stop value

Self-explanatory.

E: Run command start value is after start of run

Self-explanatory.

E: Run command stop value is before end of run

Self-explanatory.

E: Too many timesteps

The cummulative timesteps must fit in a SPARTA big integer, as as
specified by the -DSPARTA_SMALL, -DSPARTA_BIG, or -DSPARTA_BIGBIG
options in the low-level Makefile used to build SPARTA.  See Section
2.2 of the manual for details.

*/
