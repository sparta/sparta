/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(ke/molecule,ComputeKEMolecule)

#else

#ifndef DSMC_COMPUTE_KE_MOLECULE_H
#define DSMC_COMPUTE_KE_MOLECULE_H

#include "compute.h"

namespace DSMC_NS {

class ComputeKEMolecule : public Compute {
 public:
  ComputeKEMolecule(class DSMC *, int, char **);
  ~ComputeKEMolecule();
  void init();
  void compute_per_molecule();
  bigint memory_usage();

 private:
  int nmax;
  double *ke;
};

}

#endif
#endif
