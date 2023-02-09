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

#include "stdlib.h"
#include "string.h"
#include "create_box.h"
#include "domain.h"
#include "update.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

CreateBox::CreateBox(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

void CreateBox::command(int narg, char **arg)
{
  if (domain->box_exist)
    error->all(FLERR,"Cannot create_box after simulation box is defined");

  //if (domain->dimension == 2 && domain->zperiodic == 0)
  //  error->all(FLERR,"Cannot run 2d simulation with nonperiodic Z dimension");

  domain->box_exist = 1;

  if (narg != 6) error->all(FLERR,"Illegal create_box command");

  domain->boxlo[0] = atof(arg[0]);
  domain->boxhi[0] = atof(arg[1]);
  domain->boxlo[1] = atof(arg[2]);
  domain->boxhi[1] = atof(arg[3]);
  domain->boxlo[2] = atof(arg[4]);
  domain->boxhi[2] = atof(arg[5]);

  if (domain->dimension == 2) {
    if (domain->boxlo[2] >= 0.0 || domain->boxhi[2] <= 0.0)
      error->all(FLERR,
                 "Create_box z box bounds must straddle 0.0 "
                 "for 2d simulations");
  }
  if (domain->axisymmetric && domain->boxlo[1] != 0.0)
    error->all(FLERR,"Box ylo must be 0.0 for axi-symmetric model");

  // problem setup using info from header

  update->ntimestep = 0;

  domain->print_box("Created ");
  domain->set_initial_box();
  domain->set_global_box();
}
