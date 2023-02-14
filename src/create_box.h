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

#ifdef COMMAND_CLASS

CommandStyle(create_box,CreateBox)

#else

#ifndef SPARTA_CREATE_BOX_H
#define SPARTA_CREATE_BOX_H

#include "pointers.h"

namespace SPARTA_NS {

class CreateBox : protected Pointers {
 public:
  CreateBox(class SPARTA *);
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot create_box after simulation box is defined

A simulation box can only be defined once.

E: Cannot run 2d simulation with nonperiodic Z dimension

Use the boundary command to make the z dimension periodic in order to
run a 2d simulation.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Create_box z box bounds must straddle 0.0 for 2d simulations

Self-explanatory.

E: Box ylo must be 0.0 for axi-symmetric model

Self-explanatory.

*/
